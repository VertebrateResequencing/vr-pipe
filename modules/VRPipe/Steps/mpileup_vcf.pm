use VRPipe::Base;

class VRPipe::Steps::mpileup_vcf with VRPipe::StepRole {
	method options_definition {
		return { samtools_exe => VRPipe::StepOption->get(description => 'path to samtools executable',
				optional => 1,
				default_value => 'samtools'),
			   samtools_mpileup_options => VRPipe::StepOption->get(description => 'samtools mpileup options',
					   optional => 1,
					   default_value => '-C50 -aug'),
			   bcftools_exe => VRPipe::StepOption->get(description => 'path to bcftools executable',
					   optional => 1,
					   default_value => 'bcftools'),
			   bcftools_view_options => VRPipe::StepOption->get(description => 'bcftools view options',
					   optional => 1,
					   default_value => '-gcv'),
			   max_cmdline_bams => VRPipe::StepOption->get(description => 'max number of bam filenames to allow on command line', 
					   optional => 1, 
					   default_value => 10),
		};
	}
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to call variants') };
    }

    method body_sub {
        return sub {

            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            my $bcftools = $options->{bcftools_exe};
            my $bcf_view_opts = $options->{bcftools_view_options};

			my $max_cmdline_bams = $options->{max_cmdline_bams};
			if (scalar (@{$self->inputs->{bam_files}}) > $max_cmdline_bams) {
				$self->warn("[todo] Generate a bam fofn");
			}

			my $bam_list;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
				$bam_list .= "$bam_path ";
            }
			my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => "mpileup.vcf", type => 'bin');
			my $vcf_path = $vcf_file->path;

            my $req = $self->new_requirements(memory => 500, time => 1);
			my $cmd = qq[$samtools mpileup $mpileup_opts $bam_list | $bcftools view $bcf_view_opts - | bgzip -c > $vcf_path];
			$self->dispatch([$cmd, $req, {output_files => [$vcf_file]}]); 
        };
    }
    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => 'a .vcf.gz file for each set of input bam files') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run samtools mpileup for one or more bams, generating one vcf per set of bams";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
