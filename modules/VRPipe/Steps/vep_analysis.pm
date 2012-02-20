use VRPipe::Base;

class VRPipe::Steps::vep_analysis with VRPipe::StepRole {
    method options_definition {
        return { 'vep_options' => VRPipe::StepOption->get(description => 'options to vep, excluding -i or -o'),
                 'tmp_exe' => VRPipe::StepOption->get(description => 'path to script to sort on chr, temp fix for VEP sort order bug'),
                 'vep_exe' => VRPipe::StepOption->get(description => 'path to your vep executable',
                                                               optional => 1,
                                                               default_value => 'variant_effect_predictor.pl') };
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
			my $vep_exe = $options->{'vep_exe'};
			my $vep_opts = $options->{'vep_options'};
			my $tmp_exe = $options->{'tmp_exe'};
			my $cat_exe;

			if ($vep_opts =~ /-[i,o]/) {
				$self->throw("vep_options should not include the -i or -o option");
			}

			my $req = $self->new_requirements(memory => 5000, time => 1);
			foreach my $vcf_file (@{$self->inputs->{vcf_files}}) {
				my $basename = $vcf_file->basename;
				if ($basename =~ /\.vcf.gz$/) {
					$basename =~ s/\.vcf.gz$/.vep.txt/;
					$cat_exe = 'zcat';
				}
				else {
					$basename .= '.vep.txt';
					$cat_exe = 'cat';
				}
				my $vep_txt = $self->output_file(output_key => 'vep_txt', basename => $basename, type => 'txt',
													metadata => {source_vcf => $vcf_file->path->stringify});

				my $input_path = $vcf_file->path;
				my $output_path = $vep_txt->path;
				my $tmp_path = "$output_path.tmp";	# fix to vep 2.2 sort order

				my $this_cmd = "$cat_exe $input_path | $vep_exe $vep_opts -o $tmp_path; $tmp_exe $tmp_path > $output_path";

				$self->dispatch_wrapped_cmd('VRPipe::Steps::vep_analysis', 'vep_analysis', [$this_cmd, $req, {output_files => [$vep_txt]}]);
			}
		};
	}
    method outputs_definition {
        return { vep_txt => VRPipe::StepIODefinition->get(type => 'txt',
                                                             description => 'vep analysis output file',
                                                             max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Analyse VCF files with VEP";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method vep_analysis (ClassName|Object $self: Str $cmd_line) {

        my ($input_path, $output_path) = $cmd_line =~ /^\S+ (\S+) .* (\S[^;]+)$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        my $input_recs = $input_file->num_records;
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
		# Should at least be more than one output line per vcf record
        unless ($output_lines > $input_recs) {
            $output_file->unlink;
            $self->throw("VEP output has $output_lines lines, less than input vcf records $input_recs");
        }
        else {
            return 1;
        }
    }
}

1;
