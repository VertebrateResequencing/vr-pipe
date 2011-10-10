use VRPipe::Base;

class VRPipe::Steps::vcf_merge with VRPipe::StepRole {
    method options_definition {
        return { 'vcf-isec_options' => VRPipe::StepOption->get(description => '...',
							       optional => 1,
							       default_value => '-f -n +1'),
                 'vcf-isec_exe' => VRPipe::StepOption->get(description => 'path to your vcf-isec executable',
                                                           optional => 1,
                                                           default_value => 'vcf-isec'),
		 'tabix_exe' => VRPipe::StepOption->get(description => 'path to your tabix executable',
                                                        optional => 1,
                                                        default_value => 'tabix') };
    }
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'vcf',
							    description => 'compressed vcf files',
							    max_files => -1) };
    }
	method body_sub {
		return sub {
			my $self = shift;

			my $options = $self->options;
			my $tabix_exe = $options->{tabix_exe};
			my $isec_exe = $options->{'vcf-isec_exe'};
			my $isec_opts = $options->{'vcf-isec_options'};

			my $req = $self->new_requirements(memory => 500, time => 1);

			my ($merged_basename, @input_set);
			foreach my $vcf_file (@{$self->inputs->{vcf_files}}) {
				push @input_set, $vcf_file->path;

				if (!$merged_basename) {	# Use first vcf basename for output
					$merged_basename = $vcf_file->basename;
					$merged_basename =~ s/\.gz$/.merged.gz/;
                }
			}

			my $merged_vcf = $self->output_file(output_key => 'merged_vcf', basename => $merged_basename, type => 'vcf');
			my $tbi = $self->output_file(output_key => 'tbi_file', basename => "$merged_basename.tbi", type => 'bin');

			my $output_path = $merged_vcf->path;
			my $this_cmd = "$isec_exe $isec_opts @input_set | bgzip -c > $output_path; $tabix_exe -f -p vcf $output_path";

			$self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_merge', 'merge_vcf', [$this_cmd, $req, {output_files => [$merged_vcf, $tbi]}]);
		};
	}
    method outputs_definition {
        return { merged_vcf => VRPipe::StepIODefinition->get(type => 'vcf',
                                                             description => 'a merged vcf file',
                                                             max_files => 1),
                 tbi_file => VRPipe::StepIODefinition->get(type => 'bin',
                                                           description => 'a tbi file',
                                                           max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Merge soft-filtered VCFs, creating a single output VCF";
    }
    
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method merge_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($first_input_path, $output_path) = $cmd_line =~ /^vcf-isec -f -n \+1 (\S+) .* (\S[^;]+);/;
	
        my $first_input_file = VRPipe::File->get(path => $first_input_path);
        my $first_input_lines = $first_input_file->lines;
        
        $first_input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        if ($output_lines < $first_input_lines) {
            $output_file->unlink;
            $self->throw("Output VCF has $output_lines, fewer than first input $first_input_lines");
        }
        else {
            return 1;
        }
    }
}

1;
