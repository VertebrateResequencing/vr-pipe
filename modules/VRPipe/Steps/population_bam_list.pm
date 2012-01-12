use VRPipe::Base;

class VRPipe::Steps::population_bam_list with VRPipe::StepRole {
	method options_definition {
		return {
				whole_genome_mode => VRPipe::StepOption->get(description => "Indicates bam lists not split by chromosome",optional => 1,default_value => 0),
				chrom_list =>  VRPipe::StepOption->get(description => 'Names of chromosomes', optional => 1, default_value => '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y' ),
				pops_definition => VRPipe::StepOption->get(description => "populations structure definition, eg {AMR=>[qw(MXL CLM PUR)],AFR=>[qw(YRI LWK ASW)],ASN=>[qw(CHB CHS JPT)],EUR=>[qw(CEU TSI FIN GBR IBS)],}"),
		};
	}
    method inputs_definition {
        return { bam_fofn => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1, description => 'fofn of bam files, with each filename containing a population id'),
		};
    }

    method body_sub {
        return sub {

			my $self = shift;
			my $options = $self->options;
			my $whole_genome_mode = $options->{whole_genome_mode};
			my $chrom_list = $options->{chrom_list};
			my $pops_definition = $options->{pops_definition};

			my $populations;
			eval '$populations='.$pops_definition;
			$self->throw("Failed to define populations :$@") if $@;

			foreach my $bam_list (@{$self->inputs->{bam_fofn}}) {
            	my @bams;
				my $in = $bam_list->openr;
				while (my $bam=<$in>) {
					push(@bams,$bam);
				}
				close($in);

				if (!$whole_genome_mode) {
					my $chroms = [split(' ',$chrom_list)];
					for my $chr (@$chroms) {
						while (my ($pop_group,$pop_list) = each %{$populations}) {
							my $basename = "chr$chr-$pop_group";
							my $bam_fofn = $self->output_file(output_key => 'population_bam_fofn', basename => "$basename.fofn", type => 'txt');
							my $out = $bam_fofn->openw;
							foreach my $bam (@bams) {
								for my $pop (@$pop_list) {
									if ( $bam =~ /$pop/ ) {
										print $out $bam;
									}
								}
							}
							$bam_fofn->close;
						}
					}
				}
				else {
					while (my ($pop_group,$pop_list) = each %{$populations}) {
						my $basename = "pooled-$pop_group";
						my $bam_fofn = $self->output_file(output_key => 'population_bam_fofn', basename => "$basename.fofn", type => 'txt');
						my $out = $bam_fofn->openw;
						foreach my $bam (@bams) {
							for my $pop (@$pop_list) {
								if ( $bam =~ /$pop/ ) {
									print $out $bam;
								}
							}
						}
						$bam_fofn->close;
					}
				}
			}
        };
    }
    method outputs_definition {
        return { population_bam_fofn => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1, description => 'bam fofn for a population group') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Creates a seperate bam fofn for each population group, optionally by chromosome";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
