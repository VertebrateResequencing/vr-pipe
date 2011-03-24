=head1 DESCRIPTION

Just testing out DBIx::Class

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::CD extends VRPipe::Persistent {
    use VRPipe::Base::Types qw(Varchar IntSQL);
    
    has 'cdid' => (is => 'rw',
                   isa => IntSQL[16],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_auto_increment => 1,
                   is_primary_key => 1);
    
    has 'artist' => (is => 'rw',
                     isa => IntSQL[16],
                     traits => ['VRPipe::Persistent::Attributes']);
    
    has 'title' => (is => 'rw',
                    isa => Varchar[256],
                    traits => ['VRPipe::Persistent::Attributes']);
    
    has 'rank' => (is => 'rw',
                   isa => IntSQL[16],
                   traits => ['VRPipe::Persistent::Attributes']);
    
    has 'non_column_data' => (is => 'rw',
                              isa => 'Str');
    
    __PACKAGE__->make_persistent(belongs_to => {artist => 'VRPipe::Artist'},
                                 has_many => {tracks => 'VRPipe::Track'});
}

1;
