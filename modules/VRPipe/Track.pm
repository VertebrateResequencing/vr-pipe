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
                      traits => ['VRPipe::Persistent::Attributes'],
                      is_auto_increment => 1,
                      is_primary_key => 1);
    
    has 'cd' => (is => 'rw',
                 isa => IntSQL[64],
                 traits => ['VRPipe::Persistent::Attributes']);
    
    has 'title' => (is => 'rw',
                    isa => Varchar[64],
                    traits => ['VRPipe::Persistent::Attributes']);
    
    __PACKAGE__->make_persistent(belongs_to => {cd => 'VRPipe::CD'});
}

1;
