use VRPipe::Base;

class VRPipe::StepCmdSummary extends VRPipe::Persistent {
    has 'exe' => (is => 'rw',
                  isa => Varchar[255],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'version' => (is => 'rw',
                      isa => Varchar[32],
                      traits => ['VRPipe::Persistent::Attributes'],
                      is_key => 1);
    
    has 'summary' => (is => 'rw',
                      isa => Text,
                      traits => ['VRPipe::Persistent::Attributes'],
                      is_key => 1);
    
    __PACKAGE__->make_persistent();
    
    method determine_version (ClassName|Object $self: Str $cmd, Str $regex) {
        open(my $fh, "$cmd 2>&1 |") || $self->throw("Could not start $cmd");
        my $version = 0;
        
        while (<$fh>) {
            if (/$regex/) {
                $version = $1;
                last;
            }
        }
        close($fh);
        
        return $version;
    }
}

1;