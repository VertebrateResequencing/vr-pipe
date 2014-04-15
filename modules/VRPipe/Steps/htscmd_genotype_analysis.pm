
=head1 NAME

VRPipe::Steps::htscmd_genotype_analysis - a step

=head1 DESCRIPTION

The genotype check status is determined by analysing the htscmd gtcheck output
file, and computing the ratio of best-match genotype concurrence to next-best
concurrence; results are stored as metadata (gtype_analysis) on the bam file.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>, Sendu Bala <sb10@sanger.ac.uk>

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013,2014 Genome Research Limited.

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

class VRPipe::Steps::htscmd_genotype_analysis extends VRPipe::Steps::genotype_analysis {
    method description {
        return "Derive genotype check status from a htscmd gtcheck output file, computing the ratio of best-match genotype concurrence to next-best concurrence; results are stored as metadata (gtype_analysis) on the bam file.";
    }
    
    method analyse_gtcheck_output (ClassName|Object $self: Str|File :$gt_file_path!, Num :$min_ratio!, Num :$min_sites!, Num :$min_concordance!, Bool :$multiple_samples_per_individual!) {
        my $gt_file  = VRPipe::File->get(path => $gt_file_path);
        my $meta     = $gt_file->metadata;
        my $expected = $meta->{expected_sample};
        
        my $found_expected = 0;
        my ($gtype1, $score1, $gtype2, $score2, $score_expected);
        
        my $pipe = "grep -v '^#' $gt_file_path| sort -nr -k1 |"; # sort descending on concordance
        open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
        
        while (<$fh>) {
            # 0.435266        0.468085        25905095.2      25      NA20544
            my ($concordance, $uncertainty, $avg_depth, $sites, $sample) = split;
            
            next if $sites < $min_sites;
            
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
        close $fh;
        
        return $self->_analyse_output($gt_file, $min_ratio, $min_sites, $min_concordance, $multiple_samples_per_individual, $gtype1, $score1, $gtype2, $score2, $score_expected, $found_expected);
    }
}

1;
