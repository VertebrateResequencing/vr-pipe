use VRPipe::Base;

class VRPipe::Steps::vcf_concat with VRPipe::StepRole {
	method options_definition {
		return { vcf_concat_exe => VRPipe::StepOption->get(description => 'path to vcf-concat executable',
				optional => 1,
				default_value => 'vcf-concat'),
		};
	}
	method inputs_definition {
		return { vcf_files => VRPipe::StepIODefinition->get(type => 'bin', min_files => 1, max_files => -1, description => 'vcf files to concat, in chromosome order'),
		};
	}

    method body_sub {
        return sub {

            my $self = shift;
            my $options = $self->options;
            my $vcf_concat_exe = $options->{vcf_concat_exe};

			my $merge_list = $self->output_file(output_key => 'merge_list', basename => "merge_list.txt", type => 'txt', temporary => 1);
			my $ofh = $merge_list->openw;

            foreach my $vcf (@{$self->inputs->{vcf_files}}) {
				print $ofh $vcf->path,"\n";
            }

			my $concat_vcf = $self->output_file(output_key => 'concat_vcf', basename => "merged.vcf.gz", type => 'bin');
			my $concat_vcf_path = $concat_vcf->path;
			my $merge_list_path = $merge_list->path;

			my $cmd = qq[$vcf_concat_exe -f $merge_list_path | bgzip -c > $concat_vcf_path];
            my $req = $self->new_requirements(memory => 500, time => 1);
			$self->dispatch([$cmd, $req, {output_files => [$concat_vcf, $merge_list]}]); 
        };
    }
    method outputs_definition {
        return { concat_vcf => VRPipe::StepIODefinition->get(type => 'bin', max_files => 1, description => 'a merged catanated .vcf.gz file for each set of input vcfs')};
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run vcf-concat, generating one catanated vcf per input set of vcfs";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
