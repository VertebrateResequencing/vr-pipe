use VRPipe::Base;

class t::VRPipe::CD extends VRPipe::Persistent {
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
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has 'non_column_data' => (is => 'rw',
                              isa => 'Str');
    
    __PACKAGE__->make_persistent(belongs_to => {artist => 't::VRPipe::Artist'},
                                 has_many => {tracks => 't::VRPipe::Track'});
}

1;
