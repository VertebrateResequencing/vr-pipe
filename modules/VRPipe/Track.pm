=head1 DESCRIPTION

Just testing out DBIx::Class

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Track extends VRPipe::Persistent {
    use VRPipe::Base::Types qw(Varchar IntSQL);
    
    has 'trackid' => (is => 'rw',
                      isa => IntSQL[16],
                      lazy_build => 1,
                      metaclass => 'MooseX::MetaDescription::Meta::Attribute',
                      description => {
                          is_auto_increment => 1,
                          is_primary_key => 1
                      });
    
    has 'cd' => (is => 'rw',
                 lazy_build => 1,
                 isa => IntSQL[64]);
    
    has 'title' => (is => 'rw',
                    lazy_build => 1,
                    isa => Varchar[64]);
    
    __PACKAGE__->make_persistent;
    __PACKAGE__->belongs_to('cd' => 'VRPipe::CD');
}

1;
