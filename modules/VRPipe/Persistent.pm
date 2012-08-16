
=head1 NAME

VRPipe::Persistent - base class for objects that want to be persistent in the
db

=head1 SYNOPSIS
    
    use VRPipe::Base;
    
    class VRPipe::Artist extends VRPipe::Persistent {
	has 'name' => (is => 'rw',
		       isa => Varchar[64],
		       traits => ['VRPipe::Persistent::Attributes'],
		       is_key => 1);
        
	has 'agent' => (is => 'rw',
			isa => Persistent,
			coerce => 1,
			traits => ['VRPipe::Persistent::Attributes'],
			belongs_to => 'VRPipe::Agent');
        
	has 'age' => (is => 'rw',
		      isa => IntSQL[3],
		      traits => ['VRPipe::Persistent::Attributes'],
		      default => 0);
        
	has 'transient_value' => (is => 'rw', isa => 'Str');
	
	__PACKAGE__->make_persistent(has_many => [cds => 'VRPipe::CD']);
    }
    
    package main;
    
    use VRPipe::Persistent::Schema; # loads VRPipe::Artist and others
    
    # create a new artist in the db by supplying all is_keys:
    my $bob = VRPipe::Artist->create(name => 'Bob');
    
    # get an existing artist by supplying all is_keys, or just the id:
    $bob = VRPipe::Artist->get(name => 'Bob');
    $bob = VRPipe::Artist->get(id => 1);
    
    # get and simultaneously update an existing artist:
    $bob = VRPipe::Artist->get(name => 'Bob', age => 42);

get() and create() return an instance with all data columns and the full
benefit of the object, but if you need to get many objects/rows from the
database note that this is EXTREMELY SLOW. To speed up retrievals you need to
a) select all the rows of interest at once, b) retrieve only the columns of
data you're interested in, and c) avoid creation of fancy row objects and just
pull out the raw data desired.

VRPipe::Persistent is based on DBIx::Class, so multi-row selects are done using
DBIx::Class::ResultSet->search(). You can use it manually by extracting the
schema object out of a VRPipe::Persistent object, but instead it is recommended
to use one of the following convienience methods, which are all ultimately
wrappers around ResultSet->search. Where you see { ... } arguements in the
examples below, these are what you could supply as the first hashref arguement
to ResultSet->search ($cond) and are the search conditions, eg. { column_name
=> 'column_value' } to select all rows with the value 'column_value' in the
'column_name' column. \%attrs can also be supplied as the following argument,
which contain more advanced things like grouping, ordering and table joining
instructions.

If memory is not a concern (you're not getting too many rows),
get_column_values() combines the 3 speed ups in an easy-to-use method that
gives you column values you're interested in and nothing else. This is the
recommended way to do the fastest retrieval of raw data.
    
    # get_column_values() returns a list of strings if you supply a single
    # column as the first arg:
    my @names = VRPipe::Artist->get_column_values('name', { age => 30 });
    
    # or an array ref of array refs if you supply more than one column:
    my $array_ref = VRPipe::Artist->get_column_values(['name', 'agent'], { age => 30 });
    foreach my $vals (@$array_ref) {
	my ($name, $agent) = @$vals; # $agent is the id of a VRPipe::Agent 
    }
    
    # you can supply attributes (columns attribute will be overwritten) as the
    # 3rd argument:
    my ($name) = VRPipe::Artist->get_column_values('name', { age => 30 }, { rows => 1 });

If you need full instances as return values instead of raw column values, use
search(), which takes the same arguements as DBIx::Class::ResultSet->search()
and returns a list of instances in list context, but it differs by returning a
count in scalar context instead of a ResultSet:
    
    my $count = VRPipe::StepStats->search({ ... }); # very fast count(*) in SQL
    my @stepstats_instances = VRPipe::StepStats->search({ ... });
    
    # if you want a ResultSet to work with manually, use search_rs():
    my $rs = VRPipe::StepStats->search_rs({ ... }, { order_by => { -desc => ['memory'] } });
    while (my $stepstats_instance = $rs->next) {
	my $memory = $stepstats_instance->memory;
	# $rs->next loop with $stepstats_instance having all columns is slowest
    }
    # or:
    $rscolumn = $rs->get_column('memory');
    while (my $memory = $rs_column->next) {
	# $rs->next loop when only dealing with a single column is faster
    }
    # or:
    @stepstats_instances = $rs->all; # much faster than an $rs->next loop
    # or:
    my @memory_values = $rs_column->all; # fastest
    
If you're dealing with a search that could give a large number of rows and
you're worried about running out of memory, but still need greater speed than
is possible with an inefficent $rs->next loop, you can use search_paged() which
returns a VRPipe::Persistent::Pager object, which can be used like this:
    
    my $pager = VRPipe::StepStats->search_paged({ ... });
    while (my $stepstats = $pager->next) {
	# $stepstats is an array ref of VRPipe::StepStats instances
    }
    # search_paged() takes the same first two arguments as search(), but also
    # takes an optional 3rd arguement of an integer which is the number of rows
    # per page. This defaults to 5000 which stops you using unbounded memory but
    # gives almost the same efficiency and speed as search().
    
    # likewise, there is a paged version of get_column_values:
    $pager = VRPipe::StepStats->get_column_values_paged('memory', { ... });
    while (my $vals = $pager->next) {
	# $vals is an array ref of memory values
    }
    
NB: VRPipe::Persistent::Pager will reset to page 1 during a ->next loop if the
results you were searching for changes between pages. This means that what you
do during the loop should be safe to rerun on the same database row more than
once.

create() can be used to create (insert) new objects into the database, and also
to update them. However you should avoid using create() for this purpose in a
loop - it is VERY SLOW. When you don't care about the return value and just
want to insert/update rows, use bulk_create_or_update().
    
    # instead of this:
    foreach my $i (1..1000) {
	VRPipe::Job->create(cmd => "job $i", dir => '/fake_dir'); # SLOW
    }
    # collect the arguements you would have given to get() and supply them in a
    # list to bulk_create_or_update():
    my @job_args;
    foreach my $i (1..1000) {
	push(@job_args, { cmd => "job $i", dir => '/fake_dir' });
    }
    VRPipe::Job->bulk_create_or_update(@job_args); # fast

=head1 DESCRIPTION

Moose interface to DBIx::Class.

DBIx::Class::MooseColumns is OK, but I prefer my own interface here.

Use VRPipe::Base as normal to setup your class with normal 'has' sugar. For
attributes that you want stored in the database, simply specify
VRPipe::Persistent::Attributes as one of its traits. That trait will allow you
to specificy most DBIx::Class::ResultSource::add_columns args
(L<http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSource.pm#add_columns>)
as properties of your attribute. data_type is not accepted; instead your normal
'isa' determines the data_type. Your isa must be one of IntSQL, Varchar,
'Bool', Datetime, Persistent, File, Dir, FileType, AbsoluteFile, 'CodeRef',
'HashRef' or 'ArrayRef'. The last two only support simple refs with scalar
values. default_value will also be set from your attribute's default if it is
present and a simple scalar value. is_nullable defaults to false. A special
'is_key' boolean can be set which results in the column being indexed and used
as part of a multi-column (with other is_key columns) uniqueness constraint
when deciding weather to get or create a new row with get(). 'is_primary_key'
is still used to define the real key, but this is forced to be an
auto-increment called "id". 'allow_key_to_default' will allow a column to be
left out of a call to get() when that column is_key and has a default or
builder, in which case get() will behave as if you had supplied that column
with the default value for that column.

NB: for any non-persistent attributes with a default value, be sure to make
them lazy or they might not get their default values when the instances are
created.

End your class definition with a call to __PACKAGE__->make_persistent, where
you can supply the various relationship types as a hash (key as one of the
relationship methods has_many or many_to_many
(L<http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/Relationship.pm>),
value as an array ref of the args you would send to that relationship method,
or an array ref of array refs if you want to specify multiple of the same
relationship type). The other relationship types (belongs_to, has_one and
might_have) can be supplied as properties of attributes, again with an array
ref value, or just a class name string for the default configuration. You can
also supply table_name => $string if you don't want the table_name in your
database to be the same as your class basename.

For end users, create() is a convienience method that will do the equivalent of
update_or_create on a ResultSource for your class, if supplied values for all
is_key columns (with the optional exception of any allow_key_to_default
columns) and an optional instance of VRPipe::Persistent::SchemaBase to the
schema key (defaults to a production instance of VRPipe::Persistent::Schema).
You can also call get(id => $id) if you know the real auto-increment key id()
for your desired row.

clone() can be called on an instance of this class, supplying it 1 or more
is_key columns. You'll get back a (potentially) new instance with all the same
is_key column values as the instance you called clone() on, except for the
different values you supplied to clone().

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Persistent extends (DBIx::Class::Core, VRPipe::Base::Moose) { # because we're using a non-moose class, we have to specify VRPipe::Base::Moose to get Debuggable
    use MooseX::NonMoose;
    use B::Deparse;
    use Storable qw(nfreeze thaw);
    use Module::Find;
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::Persistent::ConverterFactory;
    use VRPipe::Persistent::Pager;
    use Data::Dumper;
    
    our $GLOBAL_CONNECTED_SCHEMA;
    our $deparse = B::Deparse->new("-d");
    
    our %factory_modules;
    
    __PACKAGE__->load_components(qw/InflateColumn::DateTime/);
    
    has 'id' => (
        is                => 'rw',
        isa               => IntSQL [9],
        traits            => ['VRPipe::Persistent::Attributes'],
        is_auto_increment => 1,
        is_primary_key    => 1
    );
    
    has '-result_source' => (
        is  => 'rw',
        isa => 'DBIx::Class::ResultSource::Table'
    );
    
    has '_from_non_persistent' => (
        is      => 'rw',
        isa     => 'Maybe[Object]',
        lazy    => 1,
        builder => '_determine_if_from_non_persistent'
    );
    
    # for when this instance was not retrieved via get()
    method _determine_if_from_non_persistent {
        if ($self->can('name')) {
            my $name  = $self->name;
            my $class = ref($self);
            if (exists $factory_modules{$class} && exists $factory_modules{$class}->{$name}) {
                my $factory_class = "${class}NonPersistentFactory";
                if (exists $factory_modules{factories}->{$factory_class}) {
                    eval "require $factory_class;";
                    unless ($@) {
                        return $factory_class->create($name, {});
                    }
                }
            }
        }
        
        return;
    }
    
    method make_persistent ($class: Str :$table_name?, ArrayRef :$has_many?, ArrayRef :$many_to_many?) {
        # decide on the name of the table and initialise
        unless (defined $table_name) {
            $table_name = $class;
            $table_name =~ s/.*:://;
            $table_name = lc($table_name);
        }
        $class->table($table_name);
        
        # determine what columns our table will need from the class attributes
        my (@psuedo_keys, @non_persistent, %all_attribs, %for_indexing, %indexed, %key_defaults);
        my %relationships = (belongs_to => [], has_one => [], might_have => []);
        my %flations;
        my $meta = $class->meta;
        
        my $dbtype = lc(VRPipe::Persistent::SchemaBase->get_dbtype);              # eg mysql
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {});
        
        my %ref_attribs;
        foreach my $attr ($meta->get_all_attributes) {
            my $name = $attr->name;
            $all_attribs{$name} = 1;
            unless ($attr->does('VRPipe::Persistent::Attributes')) {
                push(@non_persistent, $name);
                next;
            }
            
            my $column_info = {};
            my $vpa_meta    = VRPipe::Persistent::Attributes->meta;
            foreach my $vpa_attr ($vpa_meta->get_attribute_list) {
                next if $vpa_attr =~ /_key/;
                
                my $vpa_base = "$vpa_attr";
                $vpa_base =~ s/.*:://;
                $vpa_base = lc($vpa_base);
                if (exists $relationships{$vpa_base}) {
                    my $thing = $attr->$vpa_attr;
                    if ($thing) {
                        my $arg;
                        if (ref($thing)) {
                            $arg = $thing;
                        }
                        else {
                            $arg = [$name => $thing];
                        }
                        push(@{ $relationships{$vpa_base} }, $arg);
                    }
                    next;
                }
                
                my $predicate = $vpa_attr . '_was_set';
                next unless $attr->$predicate();
                $column_info->{$vpa_attr} = $attr->$vpa_attr;
            }
            
            if ($attr->is_primary_key) {
                $class->throw("Only id may be the primary key") unless $name eq 'id';
            }
            elsif ($attr->is_key) {
                push(@psuedo_keys, $name);
                if ($attr->allow_key_to_default) {
                    my $default = $attr->_key_default;
                    $key_defaults{$name} = ref $default ? &{$default}($class) : $default;
                }
            }
            
            # add default from our attribute if not already provided
            if (!exists $column_info->{default_value} && defined $attr->default) {
                my $default = $attr->default;
                if (ref($default)) {
                    #$default = &{$default}; #*** we assume it isn't safe to apply some dynamic value to the sql schema table definition default
                    $column_info->{is_nullable} = 1;
                }
                else {
                    $column_info->{default_value} = $default;
                }
            }
            
            # determine the type constraint that the database should use
            if ($attr->has_type_constraint) {
                my $t_c   = $attr->type_constraint;
                my $cname = $t_c->name;
                $cname =~ s/^VRPipe::.+:://;
                
                # $cname needs to be converted to something the database can use when creating the tables
                my $size       = 0;
                my $is_numeric = 0;
                if ($cname =~ /IntSQL\[(\d+)\]/) {
                    $size = $1;
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => $size, is_numeric => 1);
                }
                elsif ($cname =~ /Varchar\[(\d+)\]/) {
                    $size = $1;
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => $size, is_numeric => 0);
                }
                elsif ($cname eq 'Text') {
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => -1, is_numeric => 0);
                }
                elsif ($cname eq 'Bool') {
                    $cname = $converter->get_boolean_type();
                }
                elsif ($cname =~ /Ref/) {
                    if ($cname eq 'CodeRef') {
                        $flations{$name} = {
                            inflate => sub { eval "sub $_[0]"; },
                            deflate => sub { $deparse->coderef2text(shift); }
                        };
                    }
                    elsif ($cname eq 'HashRef[ArrayRef[Str]|HashRef[Str]|Str]|Str' || $cname eq 'ArrayRef[ArrayRef[Str]|HashRef[Str]|Str]|Str') {
                        $flations{$name} = {
                            inflate => sub { thaw(shift); },
                            deflate => sub { nfreeze(shift); }
                        };
                    }
                    elsif ($cname eq 'PersistentHashRef') {
                        $flations{$name} = {
                            inflate => sub {
                                my $hash = thaw(shift);
                                while (my ($key, $serialized) = each %$hash) { my ($class, $id) = split('~', $serialized); $hash->{$key} = $class->get(id => $id); }
                                return $hash;
                            },
                            deflate => sub {
                                my $hash = shift;
                                while (my ($key, $instance) = each %$hash) { my $ref = ref($instance); my $id = $instance->id; $hash->{$key} = "$ref~$id"; }
                                return nfreeze($hash);
                              }
                        };
                    }
                    elsif ($cname eq 'PersistentFileHashRef') {
                        $flations{$name} = {
                            inflate => sub {
                                my $hash = thaw(shift);
                                while (my ($key, $array_ref) = each %$hash) {
                                    my @array;
                                    foreach my $serialized (@$array_ref) {
                                        my ($class, $id) = split('~', $serialized);
                                        push(@array, $class->get(id => $id));
                                    }
                                    $hash->{$key} = \@array;
                                }
                                return $hash;
                            },
                            deflate => sub {
                                my $hash = shift;
                                while (my ($key, $array_ref) = each %$hash) {
                                    my @array;
                                    foreach my $instance (@$array_ref) {
                                        my $ref = ref($instance);
                                        my $id  = $instance->id;
                                        push(@array, "$ref~$id");
                                    }
                                    $hash->{$key} = \@array;
                                }
                                return nfreeze($hash);
                              }
                        };
                    }
                    elsif ($cname eq 'ArrayRefOfPersistent') {
                        $flations{$name} = {
                            inflate => sub {
                                my $array = thaw(shift);
                                my @inflated;
                                foreach my $serialized (@$array) { my ($class, $id) = split('~', $serialized); push(@inflated, $class->get(id => $id)); }
                                return \@inflated;
                            },
                            deflate => sub {
                                my $array = shift;
                                my @deflate;
                                foreach my $instance (@$array) { my $ref = ref($instance); my $id = $instance->id; push(@deflate, "$ref~$id"); }
                                return nfreeze(\@deflate);
                              }
                        };
                    }
                    else {
                        die "unsupported constraint '$cname' for attribute $name in $class\n";
                    }
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => -1, is_numeric => 0);
                    
                    $ref_attribs{$name} = 1;
                }
                elsif ($cname eq 'Datetime') {
                    $cname = $converter->get_datetime_type();
                }
                elsif ($cname eq 'Persistent') {
                    $size = 9;
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => $size, is_numeric => 1);
                    $ref_attribs{$name} = 1;
                }
                elsif ($cname eq 'File' || $cname eq 'AbsoluteFile' || $cname eq 'Dir') {
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => -1, is_numeric => 0);
                }
                elsif ($cname eq 'FileType') {
                    $size = 4;
                    ($cname, $size, $is_numeric) = $converter->get_column_info(size => $size, is_numeric => 0);
                }
                else {
                    die "unsupported constraint '$cname' for attribute $name in $class\n";
                }
                
                $column_info->{data_type}  = $cname;
                $column_info->{size}       = $size if $size;
                $column_info->{is_numeric} = $is_numeric;
            }
            else {
                die "attr $name has no constraint in $class\n";
            }
            
            # note what will need indexing by us, and seperately note everything
            # that gets indexed - DBIx::Class auto-indexes belongs_to cols
            my $bt = $attr->belongs_to;
            if ($attr->is_key) {
                $for_indexing{$name} = $column_info->{data_type} unless $bt;
                $indexed{$name} = $column_info->{data_type};
            }
            if ($bt) {
                $indexed{$name} = $column_info->{data_type};
            }
            
            # *** potential extra indexes:
            # DataElement - would query speed up with withdrawn indexed?
            # DataElementState - would query speed up with completed_steps indexed?
            # Job - do I do a blind 'running'/'finished' query?
            # PipelineSetup - worth indexing active?
            # StepState - index complete?
            # Submission - do I ever look up subs by _[sha]id? index retries, _done and _failed?
            
            # add the column in DBIx::Class, altering the name of the
            # auto-generated accessor so that we will keep our moose generated
            # accessors with their constraint checking
            my $dbic_name = '_' . $name;
            $column_info->{accessor} = $dbic_name;
            #$name =~ s/^_//; #*** no tests right now to really see if this would break things horribly; not worth risking it just to have nicer table column names
            $class->add_column($name => $column_info);
        }
        
        # set the primary key(s)
        $class->set_primary_key('id');
        
        # set relationships
        $relationships{has_many}     = $has_many     || [];
        $relationships{many_to_many} = $many_to_many || [];
        my %accessor_altered;
        foreach my $relationship (qw(belongs_to has_one might_have has_many many_to_many)) { # the order is important
            my $args = $relationships{$relationship};
            unless (ref($args->[0])) {
                $args = [$args];
            }
            
            foreach my $arg (@$args) {
                next unless @$arg;
                $class->$relationship(@$arg);
                $accessor_altered{ $arg->[0] } = 1;
            }
        }
        
        # setup inflations
        while (my ($name, $ref) = each %flations) {
            $class->inflate_column($name, $ref);
        }
        
        # now that dbic has finished creating/altering accessor methods, delete
        # them and replace with moose accessors, to give us moose type
        # constraints
        my @clear_methods;
        foreach my $attr ($meta->get_all_attributes) {
            next unless $attr->does('VRPipe::Persistent::Attributes');
            my $name         = $attr->name;
            my $dbic_name    = '_' . $name;
            my $clear_method = "_clear_$name";
            push(@clear_methods, $clear_method);
            
            if ($accessor_altered{$name}) {
                # remove DBIx::Class auto-generated accessor method
                $meta->remove_method($name);
                
                # add back the Moose accessors with their constraints etc.
                $attr->install_accessors;
            }
            
            # make the accessor get and set for DBIx::Class as well
            $meta->add_around_method_modifier(
                $name => sub {
                    my $orig = shift;
                    my $self = shift;
                    
                    # if we've been called in non-void context, we'll later want to
                    # return the current inflated dbic value if it is a ref, or the
                    # moose coerced version of it otherwise
                    my $return_val;
                    if (defined wantarray()) {
                        $return_val = exists $ref_attribs{$name} ? $self->$dbic_name() : $self->$orig();
                        unless (defined $return_val) {
                            my $dbic_val = $self->$dbic_name();
                            if (defined $dbic_val) {
                                $return_val = $self->$orig($dbic_val);
                            }
                        }
                    }
                    
                    if (@_) {
                        my $value = shift;
                        
                        if (defined $value) {
                            # first try setting in our Moose accessor, so we can see
                            # if it passes the constraint
                            $self->$orig($value);
                            
                            # now set it in the DBIC accessor
                            $self->$dbic_name($value);
                        }
                        else {
                            $self->$clear_method();
                            $self->$dbic_name(undef);
                        }
                        
                        # we deliberatly do not update in the db so that if the user
                        # is setting multiple accessors, another thread getting this
                        # object won't see a partially updated state. Users must
                        # call ->update manually (or destroy their object)
                        
                        # we do not attempt to return the set value, since we only
                        # ever return the database value, which hasn't been updated
                        # yet without that call to ->update
                    }
                    
                    return $return_val;
                }
            );
        }
        
        # like discard_changes, except we don't clumsily wipe out the whole
        # $self hash
        $meta->add_method('reselect_values_from_db' => sub { 
            my $self = shift;
            return unless $self->in_storage;
            
            if (my $current_storage = $self->get_from_storage({ force_pool => 'master' })) {
                $self->{_column_data} = $current_storage->{_column_data};
                delete $self->{_inflated_column};
                delete $self->{related_resultsets};
                bless $current_storage, 'Do::Not::Exist'; # avoid deep recursion on destruction
                
                foreach my $clear_method (@clear_methods) {
                    $self->$clear_method();
                }
            }
            else {
                $self->in_storage(0);
            }
        });
        
        # various methods need a way of getting the schema, rs, class and meta
        $meta->add_method('_class_specific' => sub { 
            my ($self, $args) = @_;
            $args ||= {};
            
            my $schema = delete $args->{schema};
            unless ($schema) {
                unless ($GLOBAL_CONNECTED_SCHEMA) {
                    eval "use VRPipe::Persistent::Schema;"; # avoid circular usage problems
                    $GLOBAL_CONNECTED_SCHEMA = VRPipe::Persistent::Schema->connect;
                }
                
                $schema = $GLOBAL_CONNECTED_SCHEMA;
            }
            
            my $rs = $schema->resultset("$class");
            
            return ($schema, $rs, $class, $meta);
        });
        
        # get() and bulk_create_or_update() need to handle search args the
        # same way
        $meta->add_method('_find_args' => sub { 
            my ($self, $args, $search_mode) = @_;
            
            my $find_args = {};
            unless ($search_mode) {
                # when doing a get(), we must supply all psuedo_keys or id
                my $id = delete $args->{id};
                if ($id) {
                    $find_args->{id} = $id;
                }
                else {
                    foreach my $key (@psuedo_keys) {
                        unless (defined $args->{$key}) {
                            if (defined $key_defaults{$key}) {
                                $find_args->{$key} = $key_defaults{$key};
                            }
                            else {
                                $self->throw("get() must be supplied all non-auto-increment keys (@psuedo_keys); missing $key");
                            }
                        }
                        else {
                            $find_args->{$key} = delete $args->{$key};
                        }
                    }
                }
            }
            
            # when dealing with refs, if we pass the ref to find(), it will
            # auto-flatten it and screw up the find call; search() will just try
            # to search for the reference and fail. find/search based on the
            # deflated string instead. Also find based on Peristent ids, not the
            # object references
            my %pks = map { $_ => 1 } @psuedo_keys;
            foreach my $hash ($search_mode ? ($args) : ($find_args, $args)) {
                while (my ($key, $val) = each %$hash) {
                    next if $key =~ /\./; # *** not sure how to handle columns from a table join
                    $self->throw("You cannot search for '$key' because $class does not have that column") unless exists $all_attribs{$key};
                    
                    unless ($search_mode) {
                        next unless exists $pks{$key}; # *** not really sure why get() doesn't work on eg. steps if we don't skip non-keys
                    }
                    
                    if ($val && ref($val)) {
                        if (UNIVERSAL::can($val, 'can')) {
                            if ($val->isa('VRPipe::Persistent')) {
                                $hash->{$key} = $val->id;
                            }
                            elsif ($val->can('stringify')) {
                                $hash->{$key} = $val->stringify;
                            }
                        }
                        elsif (exists $flations{$key}) {
                            $hash->{$key} = &{ $flations{$key}->{deflate} }($val);
                        }
                    }
                }
            }
            
            my $non_persistent_args = {};
            foreach my $key (@non_persistent) {
                exists $args->{$key} || next;
                $non_persistent_args->{$key} = delete $args->{$key};
            }
            
            return [$find_args, $non_persistent_args];
        });
        
        $meta->add_method('clone' => sub { 
            my ($self, %args) = @_;
            ref($self) || $self->throw("clone can only be called on an instance");
            $self->throw("id be supplied to clone") if $args{id};
            
            foreach my $key (@psuedo_keys) {
                unless (defined $args{$key}) {
                    $args{$key} = $self->$key();
                }
            }
            
            return $self->create(%args);
        });
        
        $meta->add_method('search_rs' => sub { 
            my ($self, $search_args, $search_attributes) = @_;
            
            my (undef, $rs) = $self->_class_specific();
            
            # when using resultset->search, a common mistake and often silently
            # dangerous bug is that you supply a value that is an instance
            # instead of an id, or you supply some other ref; _find_args will
            # fix up the search args for us
            $self->_find_args($search_args, 1);
            
            # by default we'll order by id for consistency and repeatability
            unless ($search_attributes && exists $search_attributes->{order_by}) {
                $search_attributes->{order_by} = { -asc => 'me.id' };
            }
            
            return $rs->search($search_args, $search_attributes ? $search_attributes : ());
        });
        
        $meta->add_method('_get_column_values' => sub { 
            my ($self, $column_spec, $search_args, $search_attributes) = @_;
            my @columns = ref($column_spec) ? (@$column_spec) : ($column_spec);
            $search_attributes ||= {};
            
            my $rs = $self->search_rs($search_args, { %$search_attributes, columns => \@columns });
            
            # do we have to inflate any of the columns?
            my @inflaters;
            my @is_to_inflate;
            foreach my $i (0 .. $#columns) {
                my $col = $columns[$i];
                if (exists $flations{$col}) {
                    $inflaters[$i] = $flations{$col}->{inflate};
                    push(@is_to_inflate, $i);
                }
            }
            
            my $sub;
            if (@columns == 1) {
                if (@is_to_inflate) {
                    my $inflator = $inflaters[0];
                    $sub = sub {
                        my $cursor = shift->cursor;
                        my @return;
                        foreach my $ref ($cursor->all) {
                            push(@return, &{$inflator}($ref->[0]));
                        }
                        return @return;
                    };
                }
                else {
                    $sub = sub {
                        my $cursor = shift->cursor;
                        my @return;
                        foreach my $ref ($cursor->all) {
                            push(@return, $ref->[0]);
                        }
                        return @return;
                    };
                }
            }
            else {
                if (@is_to_inflate) {
                    $sub = sub {
                        my $cursor = shift->cursor;
                        my @return;
                        foreach my $ref ($cursor->all) {
                            foreach my $i (@is_to_inflate) {
                                my $code_ref = $inflaters[$i];
                                $ref->[$i] = &{$code_ref}($ref->[$i]);
                            }
                            push(@return, $ref);
                        }
                        return \@return;
                    };
                }
                else {
                    $sub = sub {
                        return [shift->cursor->all];
                    };
                }
            }
            
            return ($rs, $sub);
        });
        
        # set up meta data to add indexes for the key columns after schema deploy
        $meta->add_attribute('cols_to_idx' => (is => 'rw', isa => 'HashRef'));
        $meta->get_attribute('cols_to_idx')->set_value($meta, \%for_indexing);
        $meta->add_attribute('idxd_cols' => (is => 'rw', isa => 'HashRef'));
        $meta->get_attribute('idxd_cols')->set_value($meta, \%indexed);
    }
    
    # txn_do with auto-retries on deadlock
    sub do_transaction {
        my ($self, $transaction, $error_message_prefix, $schema) = @_;
        $schema ||= $self->result_source->schema;
        
        my $retries     = 0;
        my $max_retries = 10;
        my $row;
        my $get_row = defined wantarray();
        while ($retries <= $max_retries) {
            try {
                if ($get_row) {
                    $row = $schema->txn_do($transaction);
                }
                else {
                    # there is a large speed increase in bulk_create_or_update
                    # if we call txn_do in void context
                    $schema->txn_do($transaction);
                }
            }
            catch ($err) {
                $self->throw("Rollback failed!") if ($err =~ /Rollback failed/);
                
                # we should attempt retries when we fail due to a create()
                # attempt when some other process got the update lock on
                # the find() first; we don't look at the $err message
                # because that may be driver-specific, and because we may
                # be nested and so $err could be unrelated
                if ($retries < $max_retries) {
                    #$self->debug("will retry get due to txn failure: $err\n");  #*** this debug call causes bizarre problems in random places
                    
                    # actually, we will check the $err message, and if it is
                    # definitely a restart situation, we won't increment
                    # retries
                    unless ($err =~ /Deadlock|try restarting transaction/) {
                        $retries++;
                    }
                    sleep(1);
                    
                    next;
                }
                else {
                    $self->throw("$error_message_prefix: $err");
                }
            }
            last;
        }
        
        if ($get_row) {
            return $row;
        }
        return;
    }
    
    # get method expects all the psuedo keys and will get or create the
    # corresponding row in the db
    sub _get {
        my ($self, $create, %args) = @_;
        
        my ($schema, $rs, $class, $meta) = $self->_class_specific(\%args);
        
        # first see if there is a corresponding $class::$name class
        my $from_non_persistent;
        if (defined $args{name} && keys %args == 1) {
            my $dir = $class . 's';
            unless (exists $factory_modules{$class}) {
                $factory_modules{$class} = { map { s/^.+:://; $_ => 1 } findallmod($dir) };
            }
            
            if (exists $factory_modules{$class}->{ $args{name} }) {
                unless (exists $factory_modules{factories}) {
                    $factory_modules{factories} = { map { $_ => 1 } grep { /NonPersistentFactory$/ } findallmod('VRPipe') };
                }
                
                my $factory_class = "${class}NonPersistentFactory";
                if (exists $factory_modules{factories}->{$factory_class}) {
                    eval "require $factory_class;";
                    die "$@\n" if $@;
                    my $obj = $factory_class->create($args{name}, {});
                    $from_non_persistent = $obj;
                    $create              = 1;
                    
                    # now setup %args based on $obj; doing things this way means
                    # we return a real Persistent object, but it is based on the
                    # very latest non-persistent code
                    my %these_args;
                    foreach my $attr ($meta->get_all_attributes) {
                        next unless $attr->does('VRPipe::Persistent::Attributes');
                        my $method = $attr->name;
                        next if $method eq 'id';
                        $these_args{$method} = $obj->$method(); # incase $method() causes the try to bomb out, we don't directly alter %args until we've finished the loop
                    }
                    while (my ($key, $val) = each %these_args) {
                        $args{$key} = $val;
                    }
                }
            }
        }
        
        my $resolve = delete $args{auto_resolve};
        if ($resolve && $self->can('resolve')) {
            my $obj = $self->_get($create, %args);
            return $obj->resolve;
        }
        
        my $fa                  = $self->_find_args(\%args);
        my %find_args           = %{ $fa->[0] || {} };
        my %non_persistent_args = %{ $fa->[1] || {} };
        
        my $transaction = sub {
            # $rs->find_or_create(...) needs all required attributes
            # supplied, but if we supply something that isn't an is_key,
            # we can generate a new row when we shouldn't do. So we
            # split up the find and create calls. Actually, search() is
            # faster than find(), and not sure we need any of the fancy
            # munging that find() does for us.
            my ($return, @extra) = $rs->search(\%find_args, { for => 'update', order_by => { -asc => 'id' } }) if keys %find_args;
            
            # there should not be any @extra, but some rare weirdness may give
            # us duplicate rows in the db; take this opportunity to delete them
            foreach my $row (@extra) {
                $row->delete;
            }
            
            if ($return) {
                # update the row with any non-key args supplied
                while (my ($method, $value) = each %args) {
                    $return->$method($value);
                }
                $return->update;
            }
            elsif ($create) {
                # create the row using all db column keys
                $return = $rs->create({ %find_args, %args });
            }
            
            if ($return) {
                # update the object with any non-persistent args supplied
                while (my ($method, $value) = each %non_persistent_args) {
                    $return->$method($value);
                }
                
                # for some reason the result_source has no schema, so
                # reattach it or inflation will break
                $return->result_source->schema($schema);
            }
            else {
                die "no matching object exists in the database";
            }
            
            return $return;
        };
        
        my @fa;
        while (my ($key, $val) = each %find_args) {
            push(@fa, "$key => $val");
        }
        my $error_message_prefix = "Failed to $class\->_get(" . join(', ', @fa) . ')';
        
        my $row = $self->do_transaction($transaction, $error_message_prefix, $schema);
        
        $row->_from_non_persistent($from_non_persistent) if $from_non_persistent;
        return $row;
    }
    
    sub get {
        my $self = shift;
        return $self->_get(0, @_);
    }
    
    sub create {
        my $self = shift;
        return $self->_get(1, @_);
    }
    
    # bulk inserts from populate() are fast, but we can easily end up with
    # duplicate rows using it, so we use a more careful but sadly slower
    # wrapper. It updates existing rows just like get().
    sub bulk_create_or_update {
        my $self = shift;
        
        my ($schema, $rs, $class) = $self->_class_specific();
        
        # when checking if we need to update instead of insert (slow), we don't
        # want to lock the table for too long, so we want to group up into
        # limited size batches
        my @find_and_nonp_args;
        my $max_per_batch = 499; # less than a second for inserts, less than 3 seconds for updates
        my $batch_num     = 0;
        my %done_args;
        my %arg_counts;
        my $arg_sets = 0;
        foreach my $args (@_) {
            # we must not have duplicate args coming through here, or we will
            # create more than 1 row with the same values, breaking our pseudo-
            # key constraints
            my $str = _args_to_str($args);
            next if exists $done_args{$str};
            $done_args{$str} = 1;
            
            my $fa = $self->_find_args($args);
            push(@{ $find_and_nonp_args[$batch_num] }, [$fa->[0], $args]);
            $batch_num++ if $#{ $find_and_nonp_args[$batch_num] } == $max_per_batch;
            
            $arg_sets++;
            while (my ($key, $val) = each %{ $fa->[0] }) {
                $arg_counts{"${key}:::$val"}++;
            }
        }
        
        my @cols_to_presearch;
        if ($arg_sets > 5) {
            while (my ($arg, $count) = each %arg_counts) {
                if ($count == $arg_sets) {
                    my ($col, $val) = split(':::', $arg);
                    push(@cols_to_presearch, ($col => $val));
                }
            }
        }
        
        my $presearch_min_id;
        foreach my $batch_i (0 .. $#find_and_nonp_args) {
            my $batch = $find_and_nonp_args[$batch_i];
            
            my $transaction = sub {
                # we do not need to lock the whole table at the start here;
                # the row-level locking we do seems to be sufficient, at least
                # for MySQL InnoDB.
                
                # we want to find existing rows, but it is slow to call
                # $rs->search on every $arg_set in @$batch; we'll first do a
                # fast check that there are any matching rows by seeing if the
                # is_key args common to all $arg_sets have a count of 0
                my $do_searches = 1;
                if (@cols_to_presearch) {
                    my @min_id_args = ();
                    if ($presearch_min_id) {
                        @min_id_args = ('me.id' => { '>' => $presearch_min_id });
                    }
                    my $count = $self->get_column_values('me.id', { @cols_to_presearch, @min_id_args }, { for => 'update', rows => $max_per_batch + 1 }); # $rs->search()->count does not lock the rows
                    if ($count == 0) {
                        $do_searches = 0;
                    }
                }
                
                # go through everything and find existing rows for update
                my %to_create;
                foreach my $arg_set (@$batch) {
                    my ($find_args, $args) = @$arg_set;
                    
                    my $return;
                    if ($do_searches) {
                        # for some reason find here is a zillion times slower than
                        # using find in get(), but using search is faster than either
                        # method in get()
                        ($return, my @extra) = $rs->search($find_args, { for => 'update' }) if keys %$find_args;
                        
                        # there should not be any @extra, but some rare
                        # weirdness may give us duplicate rows in the db; take
                        # this opportunity to delete them
                        foreach my $row (@extra) {
                            $row->delete;
                        }
                    }
                    
                    if ($return) {
                        if (keys %$args) {
                            # update the row with any non-key args supplied
                            while (my ($method, $value) = each %$args) {
                                $return->$method($value);
                            }
                            $return->update;
                        }
                    }
                    else {
                        # populate() needs everything you supply to have the
                        # same set of columns, so we must group accordingly
                        # and later call populate multiple times
                        my %populate_args = (%$find_args, %$args);
                        my $pkey = join('|', sort keys %populate_args);
                        push(@{ $to_create{$pkey} }, \%populate_args);
                    }
                }
                
                # now bulk create the rows that didn't exist
                foreach my $populate_args (values %to_create) {
                    $rs->populate($populate_args);
                }
                
                # if we have skipped checking for row existance, because our
                # transaction is for this batch only we will need to do the
                # presearch optimisation again, but obviously we just created
                # a bunch of new rows that will match the presearch, so we need
                # to get the last id of what we just created so we can ignore
                # those
                if (@cols_to_presearch && !$do_searches && $batch_i != $#find_and_nonp_args) {
                    ($presearch_min_id) = $self->get_column_values('id', {@cols_to_presearch}, { order_by => { -desc => ['id'] }, rows => 1 });
                }
            };
            
            $self->do_transaction($transaction, "Failed to $class\->bulk_create_or_update", $schema);
        }
    }
    
    sub _args_to_str {
        my $args = shift;
        my $str;
        foreach my $key (sort keys %$args) {
            my $val = $args->{$key};
            my $ref = ref($val);
            if ($ref) {
                if ($ref eq 'HASH') {
                    $val = _args_to_str($val);
                }
                elsif ($ref eq 'ARRAY') {
                    $val = join(',', sort @$val);
                }
            }
            
            $str .= "$key=>$val|";
        }
        return $str;
    }
    
    method search (ClassName|Object $self: HashRef $search_args!, HashRef $search_attributes?) {
        my $rs = $self->search_rs($search_args, $search_attributes);
        
        if (wantarray()) {
            return $rs->all;
        }
        elsif (defined wantarray()) {
            return $rs->count;
        }
        else {
            return;
        }
    }
    
    method search_paged (ClassName|Object $self: HashRef $search_args!, HashRef $search_attributes?, PositiveInt $rows_per_page?) {
        $rows_per_page ||= 5000;
        my $rs = $self->search_rs($search_args, $search_attributes);
        return VRPipe::Persistent::Pager->new(resultset => $rs, rows_per_page => $rows_per_page);
    }
    
    method get_column_values (ClassName|Object $self: Str|ArrayRef[Str] $column_spec!, HashRef $search_args!, HashRef $search_attributes?) {
        my ($rs, $sub) = $self->_get_column_values($column_spec, $search_args, $search_attributes ? $search_attributes : ());
        return &$sub($rs);
    }
    
    method get_column_values_paged (ClassName|Object $self: Str|ArrayRef[Str] $column_spec!, HashRef $search_args!, HashRef $search_attributes?, PositiveInt $rows_per_page?) {
        $rows_per_page ||= 10000;
        my ($rs, $sub) = $self->_get_column_values($column_spec, $search_args, $search_attributes ? $search_attributes : ());
        return VRPipe::Persistent::Pager->new(resultset => $rs, rows_per_page => $rows_per_page, result_method => $sub);
    }
    
    sub disconnect {
        return unless $GLOBAL_CONNECTED_SCHEMA;
        $GLOBAL_CONNECTED_SCHEMA->storage->disconnect;
    }
    
    sub reconnect {
        my $self = shift;
        return unless $GLOBAL_CONNECTED_SCHEMA;
        
        # try 10 times to connect
        my $tries     = 0;
        my $connected = 0;
        my $last_error;
        do {
            $tries++;
            eval { $GLOBAL_CONNECTED_SCHEMA->storage->ensure_connected; };
            if ($@) {
                $last_error = $@;
                warn "connection failed, $tries tries so far...\n";
                sleep(1);
            }
            else {
                $connected = 1;
            }
        } while ($tries < 10 && !$connected);
        
        unless ($connected) {
            $self->throw("$last_error");
        }
        return 1;
    }
    
    # to let you Dumper a Persistent object without tons of useless stuff
    sub _dumper_hook {
        $_[0] = bless { _strip_result_source($_[0]) }, ref($_[0]);
    }
    
    sub _strip_result_source {
        my $ref = shift;
        my %hash;
        while (my ($key, $val) = each %$ref) {
            my $val_type = ref($val);
            if ($key =~ /result_source$/) {
                $val = '* not shown for clarity *';
            }
            elsif ($val_type && ($val_type eq 'HASH' || UNIVERSAL::can($val, 'can'))) {
                $val = { _strip_result_source($val) };
            }
            $hash{$key} = $val;
        }
        return %hash;
    }
    
    method dump {
        local $Data::Dumper::Freezer = '_dumper_hook';
        Dumper($self);
    }
}

1;
