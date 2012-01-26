use VRPipe::Base;

class VRPipe::Utils::cramtools extends VRPipe::Utils::java {
    has 'cramtools_path' => (is => 'ro',
                             isa => Dir,
                             coerce => 1,
                             default => "$ENV{CRAMTOOLS}");
    
    has '+memory_multiplier' => (default => 0.9);
    
    method jar (ClassName|Object $self:) {
        return Path::Class::File->new($self->cramtools_path, 'cramtools.jar');
    }
    
    method determine_cramtools_version (ClassName|Object $self:) {
        my $cramtools_jar = $self->jar;
        my $jvm_args = $self->jvm_args(50);
        return VRPipe::StepCmdSummary->determine_version(qq[java $jvm_args -jar $cramtools_jar -h], 'v([\d\.\-]+[a-z\d]*)');
    }
}

1;
