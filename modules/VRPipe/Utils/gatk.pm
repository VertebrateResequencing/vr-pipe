use VRPipe::Base;

class VRPipe::Utils::gatk extends VRPipe::Utils::java {
    
    has 'gatk_path' => (is => 'ro',
                        isa => Dir,
                        coerce => 1,
                        default => "$ENV{GATK}");
    
    method jar (ClassName|Object $self:) {
        return Path::Class::File->new($self->gatk_path, 'GenomeAnalysisTK.jar');
    }
    method determine_gatk_version (ClassName|Object $self:) {
        my $gatk_jar = $self->jar;
        my $jvm_args = $self->jvm_args(50);
        return VRPipe::StepCmdSummary->determine_version(qq[java $jvm_args -jar $gatk_jar -h], 'v([\d\.\-]+[a-z\d]*)');
    }
}

1;
