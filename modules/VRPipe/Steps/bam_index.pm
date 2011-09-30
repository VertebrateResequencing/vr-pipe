use VRPipe::Base;

class VRPipe::Steps::bam_index with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to index') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $bai_path = $self->output_file(output_key => 'bai_files',
                                              output_dir => $bam->dir,
                                              basename => $bam->basename.'.bai',
                                              type => 'bin')->path;
                $self->dispatch([qq[$samtools index $bam_path $bai_path], $req]);
            }
        };
    }
    method outputs_definition {
        return { bai_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => 'a .bai file for each input bam file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Indexes bam files using samtools";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
