=head1 NAME

VRPipe::Persistent - base class for objects that want to be persistent in the db

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
    
    use VRPipe::Artist;
    
    # get or create a new artist in the db by supplying all is_keys:
    my $bob = VRPipe::Artist->get(name => 'Bob');

=head1 DESCRIPTION

Moose interface to DBIx::Class.

DBIx::Class::MooseColumns is OK, but I prefer my own interface here.

Use VRPipe::Base as normal to setup your class with normal 'has' sugar. For
attributes that you want stored in the database, simply specify
VRPipe::Persistent::Attributes as one of its traits. That trait will allow you
to specificy most DBIx::Class::ResultSource::add_columns args
(L<http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSource.pm#add_columns>)
as properties of your attribute. data_type is not accepted; instead your normal
'isa' determines the data_type. Your isa must be one of IntSQL, Varchar, 'Bool',
Datetime, Persistent, File, Dir, FileType, AbsoluteFile, 'CodeRef', 'HashRef' or
'ArrayRef'. The last two only support simple refs with scalar values.
default_value will also be set from your attribute's default if it is present
and a simple scalar value. is_nullable defaults to false. A special 'is_key'
boolean can be set which results in the column being indexed and used as part of
a multi-column (with other is_key columns) uniqueness constraint when deciding
weather to get or create a new row with get(). 'is_primary_key' is still used to
define the real key, but this is forced to be an auto-increment called "id".
'allow_key_to_default' will allow a column to be left out of a call to get()
when that column is_key and has a default or builder, in which case get() will
behave as if you had supplied that column with the default value for that
column.

NB: for any non-persistent attributes with a default value, be sure to make them
lazy or they might not get their default values when the instances are created.

End your class definition with a call to __PACKAGE__->make_persistent, where you
can supply the various relationship types as a hash (key as one of the
relationship methods has_many or many_to_many
(L<http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/Relationship.pm>),
value as an array ref of the args you would send to that relationship method, or
an array ref of array refs if you want to specify multiple of the same
relationship type). The other relationship types (belongs_to, has_one and
might_have) can be supplied as properties of attributes, again with an array ref
value, or just a class name string for the default configuration.
You can also supply table_name => $string if you don't want the table_name in
your database to be the same as your class basename.

For end users, get() is a convienience method that will call find_or_create on a
ResultSource for your class, if supplied values for all is_key columns (with
the optional exception of any allow_key_to_default columns) and an optional
instance of VRPipe::Persistent::SchemaBase to the schema key (defaults to a
production instance of VRPipe::Persistent::Schema). You can also call
get(id => $id) if you know the real auto-increment key id() for your desired
row.

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
    
    our $GLOBAL_CONNECTED_SCHEMA;
    our $deparse = B::Deparse->new("-d");
    
    our %factory_modules;
    
    __PACKAGE__->load_components(qw/InflateColumn::DateTime/);
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[9],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has '-result_source' => (is => 'rw',
                             isa => 'DBIx::Class::ResultSource::Table');
    
    has '_from_non_persistent' => (is => 'rw',
                                   isa => 'Maybe[Object]',
                                   lazy => 1,
                                   builder => '_determine_if_from_non_persistent');
    
    # for when this instance was not retrieved via get()
    method _determine_if_from_non_persistent {
        if ($self->can('name')) {
            my $name = $self->name;
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
        my (@psuedo_keys, @non_persistent, %for_indexing, %indexed, %key_defaults);
        my %relationships = (belongs_to => [], has_one => [], might_have => []);
        my %flations;
        my $meta = $class->meta;
	
        my $dbtype = lc(VRPipe::Persistent::SchemaBase->get_dbtype); # eg mysql
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {});
	
        foreach my $attr ($meta->get_all_attributes) {
            my $name = $attr->name;
            unless ($attr->does('VRPipe::Persistent::Attributes')) {
                push(@non_persistent, $name);
                next;
            }
            
            my $column_info = {};
            my $vpa_meta = VRPipe::Persistent::Attributes->meta;
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
                        push(@{$relationships{$vpa_base}}, $arg);
                    }
                    next;
                }
                
                my $predicate = $vpa_attr.'_was_set';
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
            if (! exists $column_info->{default_value} && defined $attr->default) {
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
                my $t_c = $attr->type_constraint;
                my $cname = $t_c->name;
                $cname =~ s/^VRPipe::.+:://;
                
                # $cname needs to be converted to something the database can use when creating the tables
                my $size = 0;
                my $is_numeric = 0;
                if ($cname =~ /IntSQL\[(\d+)\]/) {
                    $size = $1;
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> $size, is_numeric => 1);
                }
                elsif ($cname =~ /Varchar\[(\d+)\]/) {
                    $size = $1;
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> $size, is_numeric => 0);
                }
                elsif ($cname eq 'Text') {
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> -1, is_numeric => 0);
                }
                elsif ($cname eq 'Bool') {
                    $cname = $converter->get_boolean_type();
                }
                elsif ($cname =~ /Ref/) {
                    if ($cname eq 'CodeRef') {
                        $flations{$name} = { inflate => sub { eval "sub $_[0]"; },
                                             deflate => sub { $deparse->coderef2text(shift); } };
                    }
                    elsif ($cname eq 'HashRef[ArrayRef[Str]|HashRef[Str]|Str]|Str' || $cname eq 'ArrayRef[ArrayRef[Str]|HashRef[Str]|Str]|Str') {
                        $flations{$name} = { inflate => sub { thaw(shift); },
                                             deflate => sub { nfreeze(shift); } };
                    }
                    elsif ($cname eq 'PersistentHashRef') {
                        $flations{$name} = { inflate => sub { my $hash = thaw(shift); while (my ($key, $serialized) = each %$hash) { my ($class, $id) = split('~', $serialized); $hash->{$key} = $class->get(id => $id); } return $hash; },
                                             deflate => sub { my $hash = shift; while (my ($key, $instance) = each %$hash) { my $ref = ref($instance); my $id = $instance->id; $hash->{$key} = "$ref~$id"; } return nfreeze($hash); } };
                    }
                    elsif ($cname eq 'PersistentFileHashRef') {
                        $flations{$name} = { inflate => sub { my $hash = thaw(shift);
                                                              while (my ($key, $array_ref) = each %$hash) {
                                                                  my @array;
                                                                  foreach my $serialized (@$array_ref) {
                                                                      my ($class, $id) = split('~', $serialized);
                                                                      push(@array, $class->get(id => $id));
                                                                  }
                                                                  $hash->{$key} = \@array;
                                                              }
                                                              return $hash; },
                                             deflate => sub { my $hash = shift;
                                                              while (my ($key, $array_ref) = each %$hash) {
                                                                  my @array;
                                                                  foreach my $instance (@$array_ref) {
                                                                      my $ref = ref($instance);
                                                                      my $id = $instance->id;
                                                                      push(@array, "$ref~$id");
                                                                  }
                                                                  $hash->{$key} = \@array;
                                                              }
                                                              return nfreeze($hash); } };
                    }
                    elsif ($cname eq 'ArrayRefOfPersistent') {
                        $flations{$name} = { inflate => sub { my $array = thaw(shift); my @inflated; foreach my $serialized (@$array) { my ($class, $id) = split('~', $serialized); push(@inflated, $class->get(id => $id)); } return \@inflated; },
                                             deflate => sub { my $array = shift; my @deflate; foreach my $instance (@$array) { my $ref = ref($instance); my $id = $instance->id; push(@deflate, "$ref~$id"); } return nfreeze(\@deflate); } };
                    }
                    else {
                        die "unsupported constraint '$cname' for attribute $name in $class\n";
                    }
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> -1, is_numeric => 0);
                }
                elsif ($cname eq 'Datetime') {
                    $cname = $converter->get_datetime_type();
                }
                elsif ($cname eq 'Persistent') {
                    $size = 9;
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> $size, is_numeric => 1);
                }
                elsif ($cname eq 'File' || $cname eq 'AbsoluteFile' || $cname eq 'Dir') {
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> -1, is_numeric => 0);
                }
                elsif ($cname eq 'FileType') {
                    $size = 4;
                    ($cname,$size,$is_numeric) = $converter->get_column_info(size=> $size, is_numeric => 0);
                }
                else {
                    die "unsupported constraint '$cname' for attribute $name in $class\n";
                }
                
                $column_info->{data_type} = $cname;
                $column_info->{size} = $size if $size;
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
            my $dbic_name = '_'.$name;
            $column_info->{accessor} = $dbic_name;
            #$name =~ s/^_//; #*** no tests right now to really see if this would break things horribly; not worth risking it just to have nicer table column names
            $class->add_column($name => $column_info);
        }
        
        # set the primary key(s)
        $class->set_primary_key('id');
        
        # set relationships
        $relationships{has_many} = $has_many || [];
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
                $accessor_altered{$arg->[0]} = 1;
            }
        }
        
        # setup inflations
        while (my ($name, $ref) = each %flations) {
            $class->inflate_column($name, $ref);
        }
        
        # now that dbic has finished creating/altering accessor methods, delete
        # them and replace with moose accessors, to give us moose type
        # constraints
        foreach my $attr ($meta->get_all_attributes) {
            next unless $attr->does('VRPipe::Persistent::Attributes');
            my $name = $attr->name;
            my $dbic_name = '_'.$name;
            
            if ($accessor_altered{$name}) {
                # remove DBIx::Class auto-generated accessor method
                $meta->remove_method($name);
                
                # add back the Moose accessors with their constraints etc.
                $attr->install_accessors;
            }
            
            # make the accessor get and set for DBIx::Class as well
            $meta->add_around_method_modifier($name => sub {
                my $orig = shift;
                my $self = shift;
                my $moose_value = $self->$orig();
                
                # always get fresh from the db, incase another instance of this
                # row was altered and updated
                $self->reconnect;
                $self->reselect_value_from_db($name);
                
                my $dbic_value = $self->$dbic_name();
                unless (@_) {
                    # make sure we're in sync with the dbic value
                    if (defined $dbic_value) {
                        $moose_value = $self->$orig($dbic_value);
                    }
                    else {
                        $moose_value = $self->$orig();
                    }
                }
                else {
                    my $value = shift;
                    
                    if (defined $value) {
                        # first try setting in our Moose accessor, so we can see if
                        # it passes the constraint
                        $self->$orig($value);
                        
                        # now set it in the DBIC accessor
                        #$dbic_value = $self->set_column($name, $value);
                        $self->$dbic_name($value);
                    }
                    else {
                        my $clear_method = "_clear_$name";
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
                
                # if the dbic value is another Persistent object, return
                # that, otherwise prefer the moose value
                if (defined $dbic_value && ref($dbic_value)) {
                    return $dbic_value;
                }
                else {
                    return $moose_value;
                }
            });
        }
        
        # create a get method that expects all the psuedo keys and will get or
        # create the corresponding row in the db
        $meta->add_method('get' => sub {
            my ($self, %args) = @_;
            
            my $schema = delete $args{schema};
            unless ($schema) {
                unless ($GLOBAL_CONNECTED_SCHEMA) {
                    eval "use VRPipe::Persistent::Schema;"; # avoid circular usage problems
                    $GLOBAL_CONNECTED_SCHEMA = VRPipe::Persistent::Schema->connect;
                }
                
                $schema = $GLOBAL_CONNECTED_SCHEMA;
            }
            
            # first see if there is a corresponding $class::$name class
            my $from_non_persistent;
            if (defined $args{name} && keys %args == 1) {
                my $dir = $class.'s';
                unless (exists $factory_modules{$class}) {
                    $factory_modules{$class} = { map { s/^.+:://; $_ => 1 } findallmod($dir) };
                }
                
                if (exists $factory_modules{$class}->{$args{name}}) {
                    unless (exists $factory_modules{factories}) {
                        $factory_modules{factories} = { map { $_ => 1 } grep { /NonPersistentFactory$/ } findallmod('VRPipe') };
                    }
                    
                    my $factory_class = "${class}NonPersistentFactory";
                    if (exists $factory_modules{factories}->{$factory_class}) {
                        eval "require $factory_class;";
                        die "$@\n" if $@;
                        my $obj = $factory_class->create($args{name}, {});
                        $from_non_persistent = $obj;
                        
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
                my $obj = $self->get(%args);
                return $obj->resolve;
            }
            
            my %find_args;
            my $id = delete $args{id};
            if ($id) {
                %find_args = (id => $id);
            }
            else {
                foreach my $key (@psuedo_keys) {
                    unless (defined $args{$key}) {
                        if (defined $key_defaults{$key}) {
                            $find_args{$key} = $key_defaults{$key};
                        }
                        else {
                            $self->throw("get() must be supplied all non-auto-increment keys (@psuedo_keys); missing $key");
                        }
                    }
                    else {
                        $find_args{$key} = delete $args{$key};
                    }
                    
                    # when dealing with hash|array refs, if we pass the ref
                    # to find(), it will auto-flatten it and screw up the find
                    # call; find based on the frozen string instead.
                    my $val = $find_args{$key};
                    if ($val && ref($val) && (ref($val) eq 'HASH' || ref($val) eq 'ARRAY')) {
                        $find_args{$key} = nfreeze($val);
                    }
                }
            }
            
            my %non_persistent_args;
            foreach my $key (@non_persistent) {
                exists $args{$key} || next;
                $non_persistent_args{$key} = delete $args{$key};
            }
            
            $self->reconnect;
            
            my $rs = $schema->resultset("$class");
            my $row;
            
            my $transaction = sub {
                # $rs->find_or_create(...) needs all required attributes
                # supplied, but if we supply something that isn't an is_key,
                # we can generate a new row when we shouldn't do. So we
                # split up the find and create calls:
                my $return = $rs->find(\%find_args, { for => 'update' }) if keys %find_args;
                
                if ($return) {
                    # update the row with any non-key args supplied
                    while (my ($method, $value) = each %args) {
                        $return->$method($value);
                    }
                    $return->update;
                }
                else {
                    # create the row using all db column keys
                    $return = $rs->create({%find_args, %args});
                }
                
                # update the object with any non-persistent args supplied
                while (my ($method, $value) = each %non_persistent_args) {
                    $return->$method($value);
                }
                
                # for some reason the result_source has no schema, so
                # reattach it or inflation will break
                $return->result_source->schema($schema); 
                
                return $return;
            };
            
            my $retries = 0;
            my $max_retries = 10;
            while ($retries <= $max_retries) {
                try {
                    $row = $schema->txn_do($transaction);
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
                        my @fa;
                        while (my ($key, $val) = each %find_args) {
                            push(@fa, "$key => $val");
                        }
                        my $find_args = join(', ', @fa);
                        $self->throw("Failed to $class\->get($find_args): $err");
                    }
                }
                last;
            }
            
            $row->_from_non_persistent($from_non_persistent) if $from_non_persistent;
            return $row;
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
            
            return $self->get(%args);
        });
        
        # set up meta data to add indexes for the key columns after schema deploy
	$meta->add_attribute('cols_to_idx' => ( is => 'rw', isa  => 'HashRef'));
	$meta->get_attribute('cols_to_idx')->set_value($meta, \%for_indexing);
	$meta->add_attribute('idxd_cols' => ( is => 'rw', isa  => 'HashRef'));
	$meta->get_attribute('idxd_cols')->set_value($meta, \%indexed);
    }
    
    # like discard_changes, except we don't clumsily wipe out the whole $self
    # hash
    method reselect_value_from_db (Str $column) {
        return unless $self->in_storage;
        
        if (my $current_storage = $self->get_from_storage({force_pool => 'master'})) {
            $self->{_column_data}->{$column} = $current_storage->{_column_data}->{$column};
            delete $self->{_inflated_column};
            delete $self->{related_resultsets};
            bless $current_storage, 'Do::Not::Exist'; # avoid deep recursion on destruction
        }
        else {
            $self->in_storage(0);
        }
    }
    
    sub disconnect {
        return unless $GLOBAL_CONNECTED_SCHEMA;
        $GLOBAL_CONNECTED_SCHEMA->storage->disconnect;
    }
    
    sub reconnect {
        my $self = shift;
        return unless $GLOBAL_CONNECTED_SCHEMA;
        
        # try 10 times to connect
        my $tries = 0;
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
        } while ($tries < 10 && ! $connected);
        
        unless ($connected) {
            $self->throw("$last_error");
        }
        return 1;
    }
}

1;
