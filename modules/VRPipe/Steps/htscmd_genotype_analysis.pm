
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

class VRPipe::Steps::htscmd_genotype_analysis with VRPipe::StepRole {
    method options_definition {
        return {
            min_concordance => VRPipe::StepOption->create(
                description   => 'Minimum concordance to be used when applying genotype check status',
                optional      => 1,
                default_value => 0.94
            ),
            min_concordance_ratio => VRPipe::StepOption->create(
                description   => 'Minimum best-to-next-best concordance ratio to be used when applying genotype check status',
                optional      => 1,
                default_value => 1.05
            ),
            min_sites => VRPipe::StepOption->create(
                description   => 'Minimum number of sites for the data to be used in the genotype analysis',
                optional      => 1,
                default_value => 0
            ),
            multiple_samples_per_individual => VRPipe::StepOption->create(
                description   => 'If there are multiple samples per individual, set this to 1 to output the concordance before the ratio, resulting in the concordance being stored in VRTrack during auto qc, instead of the ratio',
                optional      => 1,
                default_value => 0
            )
        };
    }
    
    method inputs_definition {
        return {
            htscmd_gtcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'htscmd files containing likelihood scores for genotyping',
                metadata    => { expected_sample => 'name of expected sample', source_bam => 'input bam path' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self            = shift;
            my $options         = $self->options;
            my $min_ratio       = $options->{min_concordance_ratio};
            my $min_sites       = $options->{min_sites};
            my $min_concordance = $options->{min_concordance};
            my $mspi            = $options->{multiple_samples_per_individual};
            my $req             = $self->new_requirements(memory => 3900, time => 1);
            
            foreach my $gt_file (@{ $self->inputs->{htscmd_gtcheck_files} }) {
                my $source_bam   = $gt_file->metadata->{source_bam};
                my $gt_file_path = $gt_file->path;
                my $cmd          = "use VRPipe::Steps::htscmd_genotype_analysis; VRPipe::Steps::htscmd_genotype_analysis->analyse_htscmd_output(gt_file_path => q[$gt_file_path], source_bam => q[$source_bam], min_ratio => $min_ratio, min_sites => $min_sites, min_concordance => $min_concordance, multiple_samples_per_individual => $mspi);";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; }
    }
    
    method description {
        return "Derive genotype check status from a htscmd gtcheck output file, computing the ratio of best-match genotype concurrence to next-best concurrence; results are stored as metadata (gtype_analysis) on the bam file.";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method analyse_htscmd_output (ClassName|Object $self: Str|File :$gt_file_path!, Str|File :$source_bam!, Num :$min_ratio!, Num :$min_sites!, Num :$min_concordance!, Bool :$multiple_samples_per_individual!) {
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
                
                # if the expected sample has a score with ratio 1.000 vs score1,
                # then we'll fudge things and say gtype2 is the expected,
                # regardless of where it appeared in the file
                if (defined $gtype2 && $gtype2 ne $expected) {
                    my $ratio = $concordance != 0 ? $score1 / $concordance : $score1 / 1e-6;
                    $ratio = sprintf("%0.3f", $ratio);
                    if ($ratio == 1) {
                        $gtype2 = $sample;
                        $score2 = $concordance;
                    }
                }
            }
            
            last if ($found_expected && defined($gtype1) && defined($gtype2));
        }
        close $fh;
        
        # ratio = GT1 score / GT2 score, best-match / next-best
        my $ratio = $score2 != 0 ? $score1 / $score2 : $score1 / 1e-6;
        $ratio = sprintf("%0.3f", $ratio);
        
        # we confirm things when the top hit is much better than the second hit
        # (the ratio is higher than the min_ratio), with a special case when the
        # ratio of the scores is 1 and the concordance is over the
        # min_concordance
        
        my ($status, $exp, $con);
        if ($found_expected) {
            $exp = $expected;
            $con = $score_expected;
            
            if ($expected eq $gtype1) {
                if ($ratio >= $min_ratio || $score_expected > $min_concordance) {
                    # we're the top hit and clearly better than anything else,
                    # or have virtually no differences to the expected
                    $status = 'confirmed';
                }
                else {
                    # we could be matching any old sample within the margin of
                    # error, so call it unconfirmed
                    $status = 'unconfirmed';
                }
            }
            elsif ($expected eq $gtype2) {
                if ($ratio == 1 && $score_expected > $min_concordance) {
                    # if the scores are so close we're in the margin of error,
                    # consider the check as confirmed, but still note that the
                    # expected ($exp) and found ($gtype1) are different
                    $status = 'confirmed';
                }
                elsif ($ratio < $min_ratio) {
                    # it's pretty close to the top hit, and since it's the
                    # second hit we'll say it's not definitely wrong
                    $status = 'unconfirmed';
                }
            }
            
            # we're not in the top 2 hits, or we have a much worse score than
            # the top hit - we definitely match some unexpected sample
            $status ||= 'wrong';
        }
        else {
            $exp = 'none';
            $con = $score1;
            
            # we don't know what it's supposed to be, but ...
            if ($ratio >= $min_ratio) {
                # ... it matches clearly to the top hit and not any others
                $status = 'candidate';
            }
            else {
                # ... it's pretty similar to multiple samples, so not sure
                $status = 'unknown';
            }
        }
        
        my $gt_status = "status=$status expected=$exp found=$gtype1 " . ($multiple_samples_per_individual ? "concordance=$con ratio=$ratio" : "ratio=$ratio concordance=$con");
        print $gt_status, "\n";
        
        my $new_meta = { gtype_analysis => $gt_status };
        $gt_file->add_metadata($new_meta, replace_data => 1);
        
        my @bam_files = split(',', $source_bam);
        foreach my $bam_path (@bam_files) {
            my $bam_file = VRPipe::File->get(path => $bam_path);
            $bam_file->add_metadata($new_meta, replace_data => 1);
        }
    }
}

1;
