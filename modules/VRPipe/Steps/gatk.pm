use VRPipe::Base;

class VRPipe::Steps::gatk extends VRPipe::Steps::java {
    has 'gatk_path' => (is => 'rw',
                        isa => Dir,
                        coerce => 1);
    
    around _build_standard_options {
        return [@{$self->$orig}, 'gatk_path'];
    }
    
    our %GATK_VERSIONS;
    has 'gatk_version' => (is => 'ro',
                           isa => 'Str',
                           lazy => 1,
                           builder => 'determine_gatk_version');
    
    method jar (ClassName|Object $self:) {
        return file($self->gatk_path, 'GenomeAnalysisTK.jar');
    }
    
    method determine_gatk_version (ClassName|Object $self:) {
        my $gatk_jar = $self->jar->stringify;
        unless (defined $GATK_VERSIONS{$gatk_jar}) {
            my $jvm_args = $self->jvm_args(50);
            my $java_exe = $self->java_exe;
            $GATK_VERSIONS{$gatk_jar} = VRPipe::StepCmdSummary->determine_version(qq[$java_exe $jvm_args -jar $gatk_jar -h], 'v([\d\.\-]+[a-z\d]*)');
        }
        return $GATK_VERSIONS{$gatk_jar};
    }
    
    around options_definition {
        return { %{$self->$orig},
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
                 gatk_r_scripts => VRPipe::StepOption->get(description => 'path to GATK R script files', optional => 1, default_value => "$ENV{GATK}/resources"),
                 rscript_exe => VRPipe::StepOption->get(description => 'patht to Rscript executable', optional => 1, default_value => 'Rscript'),
                 gatk_path => VRPipe::StepOption->get(description => 'path to GATK jar files', optional => 1, default_value => "$ENV{GATK}"),
                };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub { return 1; };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generic step for using the GenomeAnalysisToolkit (GATK)";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
