=head1 DESCRIPTION

Just testing out DBIx::Class

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Artist extends VRPipe::Persistent {
    use VRPipe::Base::Types qw(Varchar IntSQL);
    
    has 'artistid' => (is => 'rw',
                       isa => IntSQL[16],
                       lazy_build => 1,
                       metaclass => 'MooseX::MetaDescription::Meta::Attribute',
                       description => {
                           is_auto_increment => 1,
                           is_primary_key => 1
                       });
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   lazy_build => 1
                   #trigger => \&_name
                   );
    
    __PACKAGE__->make_persistent;
    __PACKAGE__->has_many('cds' => 'VRPipe::CD');
    
    #__PACKAGE__->name('foo');
    #print __PACKAGE__->name, "\n";
    #before name {
    #    print "before name, name is ", $self->{name}, "\n";
    #}
    #after name {
    #    print "after name, name is ", $self->{name}, "\n";
    #}
    #method _name (Any $new_name, Any $old_name) {
    #    print "triggered name, $old_name => $new_name\n";
    #}
}

1;
