
=head1 NAME

VRPipe::Steps::bcftools_genotype_analysis - a step

=head1 DESCRIPTION

The genotype check status is determined by analysing the bcftools gtcheck
output file, and computing the ratio of best-match genotype concurrence to
next-best concurrence; results are stored as metadata (gtype_analysis) on the
bam file.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::bcftools_genotype_analysis extends VRPipe::Steps::genotype_analysis {
    method description {
        return "Derive genotype check status from a bcftools gtcheck output file, computing the ratio of best-match genotype concurrence to next-best concurrence; results are stored as metadata (gtype_analysis) on the bam file.";
    }
    
    method analyse_gtcheck_output (ClassName|Object $self: Str|File :$gt_file_path!, Num :$min_ratio!, Num :$min_sites!, Num :$min_concordance!, Bool :$multiple_samples_per_individual!) {
        my $gt_file  = VRPipe::File->get(path => $gt_file_path);
        my $meta     = $gt_file->metadata;
        my $expected = $meta->{expected_sample};
        
        my $found_expected = 0;
        my ($gtype1, $score1, $gtype2, $score2, $score_expected);
        
        my $pipe = "grep '^CN' $gt_file_path | sort -g -k3 |"; # sort ascending on discordance
        open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
        
        # we need scaled concordance values, but the output is unscaled
        # discordance, so we run through once to fix
        my ($max_disc, @data);
        while (<$fh>) {
            # CN    4.534642e+01    4.534642e+01    78  KUU25220302 0
            my (undef, undef, $discordance, $sites, $sample) = split;
            next unless $sites;
            next if $sites < $min_sites;
            
            if (!defined $max_disc || $discordance > $max_disc) {
                $max_disc = $discordance;
            }
            
            push(@data, [$discordance, $sites, $sample]);
        }
        
        if ($max_disc && @data) {
            my $min_disc     = 0;                                                                     # the minimum possible discordance is 0
            my $boundary_min = 0;
            my $boundary_max = 1 - $boundary_min;
            my $scaler       = $max_disc != $min_disc ? ($boundary_max / ($max_disc - $min_disc)) : 1;
            foreach my $datum (@data) {
                my ($discordance, $sites, $sample) = @{$datum};
                my $scaled_d = $scaler * ($discordance - $min_disc) + $boundary_min;
                my $concordance = sprintf("%0.3f", 1 - $scaled_d);
                
                if (!defined $gtype1) {
                    $gtype1 = $sample;
                    $score1 = $concordance;
                }
                elsif (!defined $gtype2) {
                    $gtype2 = $sample;
                    $score2 = $concordance;
                }
                
                if ($expected && $sample eq $expected) {
                    $found_expected = 1;
                    $score_expected = $concordance;
                    
                    # if the expected sample has a score with ratio 1.00 vs score1,
                    # then we'll fudge things and say gtype2 is the expected,
                    # regardless of where it appeared in the file
                    if (defined $gtype2 && $gtype2 ne $expected) {
                        my $ratio = $concordance != 0 ? $score1 / $concordance : $score1 / 1e-6;
                        $ratio = sprintf("%0.2f", $ratio);
                        if ($ratio == 1) {
                            $gtype2 = $sample;
                            $score2 = $concordance;
                        }
                    }
                }
                
                last if ($found_expected && defined($gtype1) && defined($gtype2));
            }
        }
        else {
            # no matches to anything; the genotype file was probably wrong, but
            # we allow this for testing purposes
            $gtype1 = 'n/a';
            $gtype2 = 'n/a';
            $score1 = '1.000';
            $score2 = '1.000';
        }
        
        close $fh;
        
        if (!$gtype2) {
            # we've only got 1 line, so the concordance will be 0; just set it
            # to 1 instead
            $score1 = '1.000';
            $score_expected = $score1 if defined $score_expected;
        }
        
        return $self->_analyse_output($gt_file, $min_ratio, $min_sites, $min_concordance, $multiple_samples_per_individual, $gtype1, $score1, $gtype2, $score2, $score_expected, $found_expected);
    }
}

1;
