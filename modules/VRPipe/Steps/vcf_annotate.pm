use VRPipe::Base;

class VRPipe::Steps::vcf_annotate with VRPipe::StepRole {
    method options_definition {
        return { 'vcf-annotate_options' => VRPipe::StepOption->get(description => 'vcf-annotate pass 1 options'),
                 'vcf-annotate_2_options' => VRPipe::StepOption->get(description => 'vcf-annotate pass 2 options'),
                 'vcf-annotate_exe' => VRPipe::StepOption->get(description => 'path to your vcf-annotate executable',
                                                               optional => 1,
                                                               default_value => 'vcf-annotate') };
    }
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'vcf',
                                                            description => 'vcf files',
                                                            max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
	    my $options = $self->options;
            my $an_exe = $options->{'vcf-annotate_exe'};
            my $an_opts = $options->{'vcf-annotate_options'};
            my $an_2_opts = $options->{'vcf-annotate_2_options'};
			my $cat_exe;
	    
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $vcf_file (@{$self->inputs->{vcf_files}}) {
                my $basename = $vcf_file->basename;
                if ($basename =~ /\.vcf.gz$/) {
                    $basename =~ s/\.vcf.gz$/.annot.vcf.gz/;
					$cat_exe = 'zcat';
                }
                else {
                    $basename =~ s/\.vcf$/.annot.vcf/;
					$cat_exe = 'cat';
                }
                my $annotated_vcf = $self->output_file(output_key => 'annotated_vcf', basename => $basename, type => 'vcf');
                
                my $input_path = $vcf_file->path;
                my $output_path = $annotated_vcf->path;
                
				# Two-pass annotation
                my $this_cmd = "$cat_exe $input_path | $an_exe $an_opts | $an_exe $an_2_opts | bgzip -c > $output_path";
		
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_annotate', 'annotate_vcf', [$this_cmd, $req, {output_files => [$annotated_vcf]}]);
            }
        };
    }
    method outputs_definition {
        return { annotated_vcf => VRPipe::StepIODefinition->get(type => 'vcf',
                                                               description => 'annotated vcf file',
                                                               max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Annotate VCF files";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method annotate_vcf (ClassName|Object $self: Str $cmd_line) {

        my ($input_path, $output_path) = $cmd_line =~ /^\S+ (\S+) .* (\S[^;]+)$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        
        my $input_lines = $input_file->lines;
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        unless ($output_lines >= $input_lines) {
            $output_file->unlink;
			$self->throw("Output VCF has fewer lines than input (input $input_lines, output $output_lines)");
        }
        else {
            return 1;
        }
    }
}

1;
