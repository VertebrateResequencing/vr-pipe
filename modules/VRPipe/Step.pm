use VRPipe::Base;

class VRPipe::Step extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'code' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    # perl -MB::Deparse -Mstrict -we 'my $deparse = B::Deparse->new("-d"); my $sub = sub { my $foo = shift; print "sub says $foo\n"; }; my $str = $deparse->coderef2text($sub); print $str, "\n"; my $ref = eval "sub $str"; &$ref("moo");'
    
    has 'description' => (is => 'rw',
                         isa => Varchar[64],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    __PACKAGE__->make_persistent();
}

1;