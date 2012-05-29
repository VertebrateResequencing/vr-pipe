use VRPipe::Base;

class VRPipe::Steps::breakdancer_bam2cfg with VRPipe::StepRole {
    method options_definition {
        return { 
            'bam2cfg_options' => VRPipe::StepOption->get(description => 'bam2cfg.pl options excluding bam file name'),
            'bam2cfg_exe' => VRPipe::StepOption->get(description => 'full path to bam2cfg.pl', optional => 1, default_value => 'bam2cfg.pl'),
        };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                            description => 'bam files',
                                                            max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            my $bam2cfg_exe = $options->{bam2cfg_exe};
            my $bam2cfg_options = $options->{'bam2cfg_options'};
            
            my $req = $self->new_requirements(memory => 500, time => 1);

            my $bam_cfg = $self->output_file(output_key => 'bam_cfg', basename => 'bd.cfg', type => 'txt');
            my $output_path = $bam_cfg->path;

			my $input_paths;
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                $input_paths .= $bam_file->path . " ";
            }
            my $cmd = "$bam2cfg_exe $bam2cfg_options $input_paths > $output_path";
            $self->dispatch_wrapped_cmd('VRPipe::Steps::breakdancer_bam2cfg', 'gen_bam_cfg', [$cmd, $req, {output_files => [$bam_cfg]}]);
        };
    }
    method outputs_definition {
        return { bam_cfg => VRPipe::StepIODefinition->get(type => 'txt',
                                                               description => 'breakdancer bam config file',
                                                               max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generates breakdancer config file from one of more input bams";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method gen_bam_cfg (ClassName|Object $self: Str $cmd_line) {

        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ /.* > (\S+)$/;
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        
        if ($output_file->num_records == 0) {
            $output_file->unlink;
			$self->throw("Output $output_path is empty)");
        }
        else {
            return 1;
        }
    }
    
}

1;
