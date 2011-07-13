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
(http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSource.pm#add_columns)
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
(http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/Relationship.pm),
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

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

class VRPipe::Persistent extends (DBIx::Class::Core, VRPipe::Base::Moose) { # because we're using a non-moose class, we have to specify VRPipe::Base::Moose to get Debuggable
    use MooseX::NonMoose;
    use B::Deparse;
    use Storable qw(nfreeze thaw);
    
    our $GLOBAL_CONNECTED_SCHEMA;
    our $deparse = B::Deparse->new("-d");
    
    __PACKAGE__->load_components(qw/InflateColumn::DateTime/);
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has '-result_source' => (is => 'rw',
                             isa => 'DBIx::Class::ResultSource::Table');
    
    method make_persistent ($class: Str :$table_name?, ArrayRef :$has_many?, ArrayRef :$many_to_many?) {
        # decide on the name of the table and initialise
        unless (defined $table_name) {
            $table_name = $class;
            $table_name =~ s/.*:://;
            $table_name = lc($table_name);
        }
        $class->table($table_name);
        
        # determine what columns our table will need from the class attributes
        my @psuedo_keys;
        my %for_indexing;
        my %key_defaults;
        my %relationships = (belongs_to => [], has_one => [], might_have => []);
        my %flations;
        my $meta = $class->meta;
        foreach my $attr ($meta->get_all_attributes) {
            next unless $attr->does('VRPipe::Persistent::Attributes');
            my $name = $attr->name;
            
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
                $for_indexing{$name} = 1;
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
                
                # $cname needs to be converted to something the database can
                # use when creating the tables, so the following cannot remain
                # hard-coded as it is now for MySQL
                my $size = 0;
                my $is_numeric = 0;
                if ($cname =~ /IntSQL\[(\d+)\]/) {
                    $cname = 'int';
                    $size = $1;
                    $is_numeric = 1;
                }
                elsif ($cname =~ /Varchar\[(\d+)\]/) {
                    $cname = 'varchar';
                    $size = $1;
                }
                elsif ($cname eq 'Bool') {
                    $cname = 'bool';
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
                    
                    $cname = 'text';
                    delete $for_indexing{$name}; # in mysql, when indexing text field we need a size, but I don't know how to specifiy the size during index creation...
                }
                elsif ($cname eq 'Datetime') {
                    $cname = 'datetime';
                }
                elsif ($cname eq 'Persistent') {
                    $cname = 'int';
                    $size = 16;
                    $is_numeric = 1;
                }
                elsif ($cname eq 'File' || $cname eq 'AbsoluteFile' || $cname eq 'Dir') {
                    $cname = 'varchar';
                    $size = 1024; #*** dangerously small?
                }
                elsif ($cname eq 'FileType') {
                    $cname = 'varchar';
                    $size = 3;
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
            
            # first see if there is a corresponding $classs::$name class
            if (defined $args{name} && keys %args == 1) {
                try {
                    eval "require ${class}s::$args{name};"; # eval within a try because can't get require to work with a variable name otherwise?!
                    die "$@\n" if $@;
                    my $factory_class = "${class}NonPersistentFactory";
                    my $obj = $factory_class->create($args{name}, {});
                    
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
                catch ($err) {
                    unless ($err =~ /^Can't locate/) {
                        $self->throw($err);
                    }
                }
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
            
            my $rs = $schema->resultset("$class");
            my $row;
            try {
                $row = $schema->txn_do(sub {
                    # $rs->find_or_create(...) needs all required attributes
                    # supplied, but if we supply something that isn't an is_key,
                    # we can generate a new row when we shouldn't do. So we
                    # split up the find and create calls:
                    my $return = $rs->find(\%find_args) if keys %find_args;
                    
                    if ($return) {
                        # update the row with any non-key args supplied
                        while (my ($method, $value) = each %args) {
                            $return->$method($value);
                        }
                        $return->update;
                    }
                    else {
                        # create the row
                        $return = $rs->create({%find_args, %args});
                    }
                    
                    # for some reason the result_source has no schema, so
                    # reattach it or inflation will break
                    $return->result_source->schema($schema); 
                    
                    return $return;
                });
            }
            catch ($err) {
                $self->throw("Rollback failed!") if ($err =~ /Rollback failed/);
                $self->throw("Failed to find_or_create: $err");
            }
            
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
        
        # add indexes for the psuedo key columns
        $meta->add_method('sqlt_deploy_hook' => sub {
            my ($self, $sqlt_table) = @_;
            $sqlt_table->add_index(name => 'psuedo_keys', fields => [keys %for_indexing]);
        });
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
    
    method disconnect {
        return unless $GLOBAL_CONNECTED_SCHEMA;
        $GLOBAL_CONNECTED_SCHEMA->storage->disconnect;
    }
}

1;
