use VRPipe::Base;

class VRPipe::Steps::vcf_index with VRPipe::StepRole {
    method options_definition {
        return {
		 'tabix_exe' => VRPipe::StepOption->get(description => 'path to tabix executable', optional => 1, default_value => 'tabix')
		};
    }
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'vcf', description => 'vcf files to index', max_files => -1) 
		};
    }
	method body_sub {
		return sub {
			my $self = shift;

			my $options = $self->options;
			my $tabix_exe = $options->{tabix_exe};

			my $req = $self->new_requirements(memory => 500, time => 1);
			foreach my $vcf (@{$self->inputs->{vcf_files}}) {

				my $vcf_path = $vcf->path;
				my $basename = $vcf->basename;
				my $tbi = $self->output_file(output_key => 'tbi_file', output_dir => $vcf->dir, basename => $basename.'.tbi', type => 'bin');

				my $cmd = "$tabix_exe -f -p vcf $vcf_path;";
				$self->dispatch([$cmd, $req, {output_files => [$tbi]}]); 
			}
		};
	}
    method outputs_definition {
        return {
                 tbi_file => VRPipe::StepIODefinition->get(type => 'bin', description => 'a tbi file', max_files => -1)
		};
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "tabix index VCF files";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
