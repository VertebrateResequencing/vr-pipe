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
    
    method make_persistent ($class: HashRef :$has_many?, HashRef :$belongs_to?) {
        # decide on the name of the table and initialise
        my $table_name = $class;
        $table_name =~ s/.*:://;
        $table_name = lc($table_name);
        $class->table($table_name);
        
        # determine what columns our table will need from the class attributes
        my @keys;
        my $meta = $class->meta;
        foreach my $attr ($meta->get_all_attributes) {
            my $name = $attr->name;
            
            my $column_info = {};
            if ($attr->does('VRPipe::Persistent::Attributes')) {
                if ($attr->is_primary_key) {
                    push(@keys, $name);
                }
                
                my $vpa_meta = VRPipe::Persistent::Attributes->meta;
                foreach my $vpa_attr ($vpa_meta->get_attribute_list) {
                    next if $vpa_attr eq 'is_primary_key';
                    my $predicate = $vpa_attr.'_was_set';
                    next unless $attr->$predicate();
                    $column_info->{$vpa_attr} = $attr->$vpa_attr;
                }
            }
            else {
                next;
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
            
            # make the accesser get and set for DBIx::Class as well
            $meta->add_around_method_modifier($name => sub {
                my $orig = shift;
                my $self = shift;
                
                #my $current_value = $self->get_column($name);
                my $current_value = $self->$dbic_name();
                return $self->$orig($current_value) unless @_;
                
                my $value = shift;
                #$self->set_column($name, $value) unless $current_value eq $value;
                $self->$dbic_name($value) unless $current_value eq $value;
                # ($self->update still needs to be called to set in db)
                
                return $self->$orig($value);
            });
        }
        
        # set the primary key(s)
        $class->set_primary_key(@keys);
        
        # set relationships
        if ($belongs_to) {
            $class->belongs_to(%{$belongs_to});
        }
        if ($has_many) {
            $class->has_many(%{$has_many});
        }
    }
}

1;
