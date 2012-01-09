use VRPipe::Base;

#Example VariantRecalibrator command - GATK v1.3
#java -Xmx4g -jar GenomeAnalysisTK.jar \
#   -T VariantRecalibrator \
#   -R reference/human_g1k_v37.fasta \
#   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
#   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
#   -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
#   -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 dbsnp_132.b37.vcf \
#   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ \
#   -recalFile path/to/output.recal \
#   -tranchesFile path/to/output.tranches \
#   -rscriptFile path/to/output.plots.R

class VRPipe::Steps::gatk_recalibrate_variants extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{$self->$orig},
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to reference genome fasta'),
                 var_recal_opts => VRPipe::StepOption->get(description => 'options for GATK VariantRecalibrator, excluding reference genome, input and output'),
               };
    }

    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => 'one or more tabixed vcf files for variant recalibration'),
		};
    }

    method body_sub {
        return sub {
            use VRPipe::Utils::gatk;
            
            my $self = shift;
            my $options = $self->options;
            my $gatk = VRPipe::Utils::gatk->new(gatk_path => $options->{gatk_path}, java_exe => $options->{java_exe});
            
			my $reference_fasta = $options->{reference_fasta};
            my $var_recal_opts = $options->{var_recal_opts};

            my $req = $self->new_requirements(memory => 1200, time => 1);
            my $jvm_args = $gatk->jvm_args($req->memory);

            foreach my $vcf (@{$self->inputs->{vcf_files}}) {
                my $vcf_path = $vcf->path;
				my $basename = $vcf->basename;
				if ($basename =~ /\.vcf.gz$/) {
					$basename =~ s/\.vcf.gz$/.recal.csv/;
				}
				else {
					$basename =~ s/\.vcf$/.recal.csv/;
				}

				my $recal_file = $self->output_file(output_key => 'recal_files', basename => $basename, type => 'txt');
				my $recal_path = $recal_file->path;
				$basename =~ s/recal\.csv/tranches/;
				my $tranches_file = $self->output_file(output_key => 'tranches_files', basename => $basename, type => 'txt');
				my $tranches_path = $tranches_file->path;

				my $cmd = $gatk->java_exe.qq[ $jvm_args -jar ].$gatk->jar.qq[ -T VariantRecalibrator -R $reference_fasta -input $vcf_path -recalFile $recal_path -tranchesFile $tranches_path $var_recal_opts ];
#				$self->warn($cmd);
				$self->dispatch([$cmd, $req, {output_files => [$recal_file, $tranches_file]}]); 
			}
        };
    }
    method outputs_definition {
        return {
			recal_files => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1, description => 'a recalibration table file in CSV format for each input vcf'),
			tranches_files => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1, description => 'a tranches file for each vcf used by ApplyRecalibration'),
		};
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run gatk VariantRecalibrator on vcf files, generating recalibration files for use by ApplyRecalibration walker";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
