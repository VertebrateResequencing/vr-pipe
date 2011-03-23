=head1 DESCRIPTION

Moose interface to DBIx::Class

DBIx::Class::MooseColumns is OK, but I prefer my own interface here.


=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Persistent extends DBIx::Class::Core {
    use MooseX::NonMoose;
    use VRPipe::Base::Types qw(Varchar IntSQL);
    
    method make_persistent ($class:) {
        # decide on the name of the table and initialise
        my $table_name = $class;
        $table_name =~ s/.*:://;
        $table_name = lc($table_name);
        $class->table($table_name);
        
        # determine what columns our table will need from the class attributes
        my @keys;
        my $meta = $class->meta;
        for my $attr ($meta->get_all_attributes) {
            my $name = $attr->name;
            
            my $column_info = {};
            if ($attr->can('description')) {
                my $desc = $attr->description;
                if ($desc) {
                    # we have 2 special keys is_primary_key and is_transient
                    my $is_transient = delete $desc->{is_transient};
                    if ($is_transient) {
                        next;
                    }
                    
                    my $is_key = delete $desc->{is_primary_key};
                    push(@keys, $name) if $is_key;
                    
                    while (my ($key, $value) = each %{$desc}) {
                        $column_info->{$key} = $value;
                    }
                }
            }
            
            # we default to is_nullable => 0
            unless (defined $column_info->{is_nullable}) {
                $column_info->{is_nullable} = 0;
            }
            
            # set the type constraint
            if ($attr->has_type_constraint) {
                my $t_c = $attr->type_constraint;
                my $cname = $t_c->name;
                
                # $cname needs to be converted to something the database can
                # use when creating the tables, so the following cannot remain
                # hard-coded as it is now for MySQL
                my $size = 0;
                if ($cname =~ /IntSQL\[(\d+)\]/) {
                    $cname = 'int';
                    $size = $1;
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
            }
            else {
                die "attr $name has no constraint in $class\n";
            }
            
            # add the column in DBIx::Class
            $class->add_column($name => $column_info);
            
            # remove DBIx::Class auto-generated accessor method
            $meta->remove_method($name);
            
            # add back the Moose accessors with their constraints etc.
            $attr->install_accessors;
            
            # make the accesser get and set for DBIx::Class as well
            $meta->add_around_method_modifier($name => sub {
                my $orig = shift;
                my $self = shift;
                
                my $current_value = $self->get_column($name);
                return $self->$orig($current_value) unless @_;
                
                my $value = shift;
                $self->set_column($name, $value) unless $current_value eq $value;
                # ($self->update still needs to be called to set in db)
                
                return $self->$orig($value);
            });
        }
        
        # set the primary key(s)
        $class->set_primary_key(@keys);
    }
}

1;
