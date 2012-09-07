use VRPipe::Base;

class VRPipe::Steps::conifer_analyze with VRPipe::StepRole {
    method options_definition {
        return {
            python_exe  => VRPipe::StepOption->create(description => 'full path to python executable',               optional => 1, default_value => 'python'),
            conifer_py  => VRPipe::StepOption->create(description => 'full path to conifer.py',                      optional => 1, default_value => 'conifer.py'),
            probes_file => VRPipe::StepOption->create(description => 'probes / target definition file'),
            write_svals => VRPipe::StepOption->create(description => 'output singular values from the SVD analysis', optional => 1, default_value => 1),
            write_sd    => VRPipe::StepOption->create(description => 'output standard deviation for each sample',    optional => 1, default_value => 1),
            svd_remove  => VRPipe::StepOption->create(description => 'the number of SVD components to remove',       optional => 1, default_value => 0),
        };
    }
    
    method inputs_definition {
        return {
            rpkm_in => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'set of conifer rpkm text files to be analyzed as a group',
                max_files   => -1
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options     = $self->options;
            my $python_exe  = $options->{python_exe};
            my $conifer_py  = $options->{conifer_py};
            my $probes_file = $options->{probes_file};
            my $write_svals = $options->{'write_svals'};
            my $write_sd    = $options->{'write_sd'};
            my $svd_remove  = $options->{'svd_remove'};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my $rpkm_file = $self->inputs->{rpkm_in}[0]; # get rpkm dir from first input file
            my $rpkm_dir  = $rpkm_file->dir;
            
            my $analysis_hdf5 = $self->output_file(output_key => 'analysis_hdf5', basename => 'analysis.hdf5', type => 'bin');
            my $output_path = $analysis_hdf5->path;
            
            my $cmd = "$python_exe $conifer_py analyze --probes $probes_file --rpkm_dir $rpkm_dir --output $output_path";
            $cmd .= " --svd $svd_remove";
            $cmd .= " --write_svals " . $analysis_hdf5->dir . "/singular_values.txt" if $write_svals;
            $cmd .= " --write_sd " . $analysis_hdf5->dir . "/sd_values.txt"          if $write_sd;
            $self->warn($cmd);
            
            $self->dispatch_wrapped_cmd('VRPipe::Steps::conifer_analyze', 'run_analysis', [$cmd, $req, { output_files => [$analysis_hdf5] }]);
        };
    }
    
    method outputs_definition {
        return {
            analysis_hdf5 => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'HDF5 file of SVD-ZRPKM values for each sample in rpkm set',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs conifer analysis for a set of rpkm input files, generating HDF5 file of SVD-ZRPKM values for each sample in the set";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method run_analysis (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ / --output (\S+)/;
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        
        if ($output_file->lines == 0) {
            $output_file->unlink;
            $self->throw("Output $output_path is empty)");
        }
        else {
            return 1;
        }
    }

}

1;
