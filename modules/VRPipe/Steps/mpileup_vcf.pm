use VRPipe::Base;

class VRPipe::Steps::mpileup_vcf extends VRPipe::Steps::mpileup_bcf {

    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            my $reference_fasta = $options->{reference_fasta};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            my $interval_list = $options->{interval_list};
            my $bcftools = $options->{bcftools_exe};
            my $bcf_view_opts = $options->{bcftools_view_options};

			my $max_cmdline_bams = $options->{max_cmdline_bams};
			if (scalar (@{$self->inputs->{bam_files}}) > $max_cmdline_bams) {
				$self->warn("[todo] Generate a bam fofn");
			}

            my $req = $self->new_requirements(memory => 500, time => 1);
			my $bam_list;
			my ($bam_metadata,$basename);

			# if more than one bam, vcf basename and any meta data will be based upon the last one
            foreach my $bam (@{$self->inputs->{bam_files}}) {
				$bam_metadata = $bam->metadata;	
				$basename = $bam->basename;	
                my $bam_path = $bam->path;
				$bam_list .= "$bam_path ";
            }

			$mpileup_opts .= " -l $interval_list " if $interval_list;
			$basename =~ s/\.bam/.mpileup.vcf.gz/;

			my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => $basename, type => 'vcf');
			my $vcf_path = $vcf_file->path;
			if ($bam_metadata) {
				$vcf_file->add_metadata($bam_metadata);
			}

			my $cmd = qq[$samtools mpileup $mpileup_opts -f $reference_fasta $bam_list | $bcftools view $bcf_view_opts - | bgzip -c > $vcf_path];
			$self->dispatch([$cmd, $req, {output_files => [$vcf_file]}]); 
	    };
    }
    method description {
        return "Run samtools mpileup and bcftools for one or more bams, generating one vcf without an intermediate bcf";
    }

    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'vcf', max_files => -1, description => 'a vcf file for each set of one or more input bams') };
    }

}

1;
