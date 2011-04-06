=head1 NAME

VRPipe::Persistent - base class for objects that want to be persistent in the db

=head1 SYNOPSIS

use VRPipe::Base;

class VRPipe::Artist extends VRPipe::Persistent {
    has 'artistid' => (is => 'rw',
                       isa => IntSQL[16],
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_auto_increment => 1,
                       is_primary_key => 1);

    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_primary_key => 1);

    has 'age' => (is => 'rw',
                  isa => IntSQL[3],
                  traits => ['VRPipe::Persistent::Attributes'],
                  default => 0);

    has 'transient_value' => (is => 'rw', isa => 'Str');
    
    __PACKAGE__->make_persistent(has_many => [cds => 'VRPipe::CD']);
}

package main;

use VRPipe::Artist;

# get or create a new artist in the db by supplying at least all primary keys
# (except auto_increment keys):
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
'isa' determines the data_type. Your isa must be one of IntSQL|Varchar.
default_value will also be set from your attribute's default if it is present
and a simple scalar value. is_nullable defaults to false.

End your class definition with a call to __PACKAGE__->make_persistent, where you
can supply the various relationship types as a hash (key as one of the
relationship methods
(http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/Relationship.pm),
value as an array ref of the args you would send to that relationship method),
and also table_name => $string if you don't want the table_name in your database
to be the same as your class basename.

For end users, get() is a convience method that will call find_or_create on a
ResultSource for your class, if supplied values for at least all key columns
(except auto_increment columns) and an instance of VRPipe::Persistent::Schema to
the schema key.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

class VRPipe::Persistent extends (DBIx::Class::Core, VRPipe::Base::Moose) { # because we're using a non-moose class, we have to specify VRPipe::Base::Moose to get Debuggable
    use MooseX::NonMoose;
    use VRPipe::Persistent::Schema;
    
    has '-result_source' => (is => 'rw', isa => 'DBIx::Class::ResultSource::Table');
    
    method make_persistent ($class: Str :$table_name?, ArrayRef :$has_many?, ArrayRef :$has_one?, ArrayRef :$belongs_to?, ArrayRef :$might_have?, ArrayRef :$many_to_many?) {
        # decide on the name of the table and initialise
        unless (defined $table_name) {
            $table_name = $class;
            $table_name =~ s/.*:://;
            $table_name = lc($table_name);
        }
        $class->table($table_name);
        
        # determine what columns our table will need from the class attributes
        my @keys;
        my @non_auto_keys;
        my $meta = $class->meta;
        foreach my $attr ($meta->get_all_attributes) {
            my $name = $attr->name;
            
            my $column_info = {};
            if ($attr->does('VRPipe::Persistent::Attributes')) {
                my $vpa_meta = VRPipe::Persistent::Attributes->meta;
                foreach my $vpa_attr ($vpa_meta->get_attribute_list) {
                    next if $vpa_attr eq 'is_primary_key';
                    my $predicate = $vpa_attr.'_was_set';
                    next unless $attr->$predicate();
                    $column_info->{$vpa_attr} = $attr->$vpa_attr;
                }
                
                if ($attr->is_primary_key) {
                    push(@keys, $name);
                    push(@non_auto_keys, $name) unless (exists $column_info->{is_auto_increment} && $column_info->{is_auto_increment});
                }
            }
            else {
                next;
            }
            
            # add default from our attribute if not already provided
            if (! exists $column_info->{default_value} && defined $attr->default) {
                $column_info->{default_value} = $attr->default;
            }
            
            # determine the type constraint that the database should use
            if ($attr->has_type_constraint) {
                my $t_c = $attr->type_constraint;
                my $cname = $t_c->name;
                
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
            # accessors with their contraint checking
            my $dbic_name = '_'.$name;
            $column_info->{accessor} = $dbic_name;
            $class->add_column($name => $column_info);
            
            # remove DBIx::Class auto-generated accessor method
            #$meta->remove_method($name);
            
            # add back the Moose accessors with their constraints etc.
            #$attr->install_accessors;
            
            # make the accessor get and set for DBIx::Class as well
            $meta->add_around_method_modifier($name => sub {
                my $orig = shift;
                my $self = shift;
                
                #my $current_value = $self->get_column($name);
                my $current_value = $self->$dbic_name();
                unless (@_) {
                    return defined $current_value ? $self->$orig($current_value) : $self->$orig();
                }
                
                my $value = shift;
                #$self->set_column($name, $value) unless $current_value eq $value;
                unless ($current_value eq $value) {
                    # first try setting in our Moose accessor, so we can see if
                    # it passes the constraint
                    my $return = $self->$orig($value);
                    
                    # now set it in the DBIC accessor
                    $self->$dbic_name($value);
                    
                    # we deliberatly do not update in the db so that if the user
                    # is setting multiple accessors, another thread getting this
                    # object won't see a partially updated state. Users must
                    # call ->update manually (or destroy their object)
                    #$self->update;
                    
                    return $return;
                }
                
                return $value;
            });
        }
        
        # set the primary key(s)
        $class->set_primary_key(@keys);
        
        # set relationships
        if ($belongs_to) {
            $class->belongs_to(@{$belongs_to});
        }
        if ($has_many) {
            $class->has_many(@{$has_many});
        }
        if ($has_one) {
            $class->has_one(@{$has_one});
        }
        if ($might_have) {
            $class->might_have(@{$might_have});
        }
        if ($many_to_many) {
            $class->many_to_many(@{$many_to_many});
        }
        
        # create a get method that expects all the primary keys and will get or
        # create the corresponding row in the db
        $meta->add_method('get' => sub {
            my ($self, %args) = @_;
            my $schema = delete $args{schema} || VRPipe::Persistent::Schema->connect;
            foreach my $key (@non_auto_keys) {
                unless (defined $args{$key}) {
                    $self->throw("get() must be supplied all non-auto-increment keys (@non_auto_keys); missing $key");
                }
            }
            
            my $rs = $schema->resultset("$class");
            return $rs->find_or_create(%args);
        });
    }
    
    method DEMOLISH {
        $self->update;
    }
}

1;
