use VRPipe::Base;

class VRPipe::Steps::mpileup_vcf with VRPipe::StepRole {
	method options_definition {
		return { samtools_exe => VRPipe::StepOption->get(description => 'path to samtools executable',
				optional => 1,
				default_value => 'samtools'),
			   samtools_mpileup_options => VRPipe::StepOption->get(description => 'samtools mpileup options excluding -f',
					   optional => 1,
					   default_value => '-C50 -aug'),
			   reference_fasta => VRPipe::StepOption->get(description => 'absolute path to reference genome fasta', 
					   optional => 1, 
					   default_value => '/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.fasta'),
			   bcftools_exe => VRPipe::StepOption->get(description => 'path to bcftools executable',
					   optional => 1,
					   default_value => 'bcftools'),
			   bcftools_view_options => VRPipe::StepOption->get(description => 'bcftools view options',
					   optional => 1,
					   default_value => '-gcv'),
			   max_cmdline_bams => VRPipe::StepOption->get(description => 'max number of bam filenames to allow on command line', 
					   optional => 1, 
					   default_value => 10),
			   interval_list => VRPipe::StepOption->get(description => 'absolute path to targets interval list file for -l option', 
					   optional => 1,),
		};
	}
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to call variants'),
		chunked_regions_file => VRPipe::StepIODefinition->get(type => 'txt', min_files => 0, max_files => 1, description => 'file of chomosome region chunks to run concurrently'),
		 };
    }

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
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
				$bam_list .= "$bam_path ";
            }

			# perform concurrent analyses if optional chunk file is present
			if  ($self->inputs->{chunked_regions_file}) {
				my $chunk_file = $self->inputs->{chunked_regions_file}[0];
				my $cfh = $chunk_file->openr;
				while (<$cfh>) {
					my ($chr,$from,$to) = split;
					my $region = "${chr}:${from}-${to}";
					my $chunk_opts = "$mpileup_opts -r ${chr}:${from}-${to}";

					my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => "mpileup_par_${chr}_${from}_${to}.vcf.gz", type => 'bin');
					my $vcf_path = $vcf_file->path;

					my $cmd = qq[$samtools mpileup $mpileup_opts -f $reference_fasta $chunk_opts $bam_list | $bcftools view $bcf_view_opts - | bgzip -c > $vcf_path];
					$self->dispatch([$cmd, $req, {output_files => [$vcf_file]}]); 
				}
			}
			else {
				$mpileup_opts .= " -l $interval_list " if $interval_list;

				my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => "mpileup.vcf.gz", type => 'bin');
				my $vcf_path = $vcf_file->path;

				my $cmd = qq[$samtools mpileup $mpileup_opts -f $reference_fasta $bam_list | $bcftools view $bcf_view_opts - | bgzip -c > $vcf_path];
				$self->dispatch([$cmd, $req, {output_files => [$vcf_file]}]); 
			}
        };
    }
    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => 'either a single .vcf.gz file, or a chunk set of vcf.gz files, for each set of input bam files') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run samtools mpileup for one or more bams, generating either one vcf per set of bams, or a chunk set per bam set if chunked_regions_file is provided";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
