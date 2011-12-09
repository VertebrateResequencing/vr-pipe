use VRPipe::Base;

#Example generic command for multi-sample SNP calling
# java -jar GenomeAnalysisTK.jar \
#   -R resources/Homo_sapiens_assembly18.fasta \
#   -T UnifiedGenotyper \
#   -I sample1.bam [-I sample2.bam ...] \
#   -o snps.raw.vcf \
#   -stand_call_conf [50.0] \
#   -stand_emit_conf 10.0 \
#   -dcov [50] \

class VRPipe::Steps::gatk_genotype extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{$self->$orig},
                 genotyper_opts => VRPipe::StepOption->get(description => 'options for GATK UnifiedGenotyper, excluding -R,-D,-I,-o'),
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to reference genome fasta', 
					optional => 1, 
					default_value => '/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.fasta'),
                 max_cmdline_bams => VRPipe::StepOption->get(description => 'max number of bam filenames to allow on command line', 
					optional => 1, 
					default_value => 10),
                 interval_list => VRPipe::StepOption->get(description => 'absolute path to targets interval list file for -L option', 
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
            use VRPipe::Utils::gatk;
            
            my $self = shift;
            my $options = $self->options;
            my $gatk = VRPipe::Utils::gatk->new(gatk_path => $options->{gatk_path}, java_exe => $options->{java_exe});
            
            my $reference_fasta = $options->{reference_fasta};
            my $genotyper_opts = $options->{genotyper_opts};
            my $max_cmdline_bams = $options->{max_cmdline_bams};

            my $interval_list = $options->{interval_list};
            
            my $req = $self->new_requirements(memory => 1200, time => 1);
            my $jvm_args = $gatk->jvm_args($req->memory);

			if (scalar (@{$self->inputs->{bam_files}}) > $max_cmdline_bams) {
				$self->warn("[todo] Generate a bam fofn");
			}
            
			my $bam_list;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
				$bam_list .= "-I $bam_path ";
            }
			
			# perform concurrent analyses if optional chunk file is present
			if  ($self->inputs->{chunked_regions_file}) {
				my $chunk_file = $self->inputs->{chunked_regions_file}[0];
				my $cfh = $chunk_file->openr;
				while (<$cfh>) {
					my ($chr,$from,$to) = split;
					my $region = "${chr}:${from}-${to}";
					my $chunk_opts = "$genotyper_opts -L ${chr}:${from}-${to}";

					my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => "gatk_par_${chr}_${from}_${to}.vcf.gz", type => 'bin');
					my $vcf_path = $vcf_file->path;

					my $cmd = $gatk->java_exe.qq[ $jvm_args -jar ].$gatk->jar.qq[ -T UnifiedGenotyper -R $reference_fasta $chunk_opts $bam_list -o $vcf_path ];
					$self->dispatch([$cmd, $req, {output_files => [$vcf_file]}]); 
				}
			}
			else {
				$genotyper_opts .= " -L $interval_list " if $interval_list;
				my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => "gatk_var.vcf.gz", type => 'bin');
				my $vcf_path = $vcf_file->path;

				$self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'GenomeAnalysisTK', 
					version => $gatk->determine_gatk_version(),
					summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference_fasta $genotyper_opts -I $bam_path -o $vcf_path'));

				my $cmd = $gatk->java_exe.qq[ $jvm_args -jar ].$gatk->jar.qq[ -T UnifiedGenotyper -R $reference_fasta $genotyper_opts $bam_list -o $vcf_path ];
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
        return "Run gatk UnifiedGenotyper for one or more bams, generating either one vcf per set of bams, or a chunk set per bam set if chunked_regions_file is provided";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
