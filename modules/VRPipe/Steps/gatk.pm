use VRPipe::Base;

class VRPipe::Steps::gatk extends VRPipe::Steps::java {
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
