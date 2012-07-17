=head1 NAME

VRPipe::Steps::population_bam_list - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Steps::population_bam_list with VRPipe::StepRole {
	method options_definition {
		return {
				whole_genome_mode => VRPipe::StepOption->create(description => "Indicates bam lists not split by chromosome",optional => 1,default_value => 0),
				chrom_list =>  VRPipe::StepOption->create(description => 'Names of chromosomes', optional => 1, default_value => '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y' ),
				pops_definition => VRPipe::StepOption->create(description => "populations structure definition, eg {AMR=>[qw(MXL CLM PUR)],AFR=>[qw(YRI LWK ASW)],ASN=>[qw(CHB CHS JPT)],EUR=>[qw(CEU TSI FIN GBR IBS)],}"),
		};
	}
    method inputs_definition {
        return { bam_fofn => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'fofn of bam files, with each filename containing a population id'),
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
        return { population_bam_fofn => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'bam fofn for a population group') };
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
