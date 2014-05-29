
=head1 NAME

VRPipe::Steps::genotype_analysis - a step

=head1 DESCRIPTION

The genotype check status is determined by analysing the gtcheck output file,
and computing the ratio of best-match genotype concurrence to next-best
concurrence; results are stored as metadata (gtype_analysis) on the bam file.

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

class VRPipe::Steps::genotype_analysis with VRPipe::StepRole {
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
            gtcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'gtcheck files containing likelihood scores for genotyping',
                metadata    => { expected_sample => 'name of expected sample', source_bam => 'input bam path' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self            = shift;
            my $class           = ref($self);
            my $options         = $self->options;
            my $min_ratio       = $options->{min_concordance_ratio};
            my $min_sites       = $options->{min_sites};
            my $min_concordance = $options->{min_concordance};
            my $mspi            = $options->{multiple_samples_per_individual};
            my $req             = $self->new_requirements(memory => 3900, time => 1);
            
            foreach my $gt_file (@{ $self->inputs->{gtcheck_files} }) {
                my $source_bam   = $gt_file->metadata->{source_bam};
                my $gt_file_path = $gt_file->path;
                my $cmd          = "use $class; $class->analyse_gtcheck_output(gt_file_path => q[$gt_file_path], min_ratio => $min_ratio, min_sites => $min_sites, min_concordance => $min_concordance, multiple_samples_per_individual => $mspi);";
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
        return "Generic step for other genotype analysis steps to inherit from.";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method _analyse_output (ClassName|Object $self: VRPipe::File $gt_file!, Num $min_ratio!, Num $min_sites!, Num $min_concordance!, Bool $multiple_samples_per_individual!, Str $gtype1!, Num $score1!, Maybe[Str] $gtype2?, Maybe[Num] $score2?, Maybe[Num] $score_expected?, Bool $found_expected?) {
        my $meta       = $gt_file->metadata;
        my $expected   = $meta->{expected_sample};
        my $source_bam = $meta->{source_bam};
        
        my ($ratio, $status, $exp, $con);
        if ($gtype2) {
            # ratio = GT1 score / GT2 score, best-match / next-best
            $ratio = $score2 != 0 ? $score1 / $score2 : $score1 / 1e-6;
            $ratio = sprintf("%0.2f", $ratio);
            
            # we confirm things when the top hit is much better than the second hit
            # (the ratio is higher than the min_ratio), with a special case when the
            # ratio of the scores is 1 and the concordance is over the
            # min_concordance
            
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
        }
        else {
            # this edge case probably only applies in testing: the gtcheck file
            # only had 1 row with enough sites
            $ratio = "1.00";
            if ($found_expected) {
                $status = 'confirmed';
                $exp    = $expected;
                $con    = $score_expected;
            }
            else {
                $status = 'unknown';
                $exp    = 'none';
                $con    = $score1;
            }
        }
        
        my $gt_status = "status=$status expected=$exp found=$gtype1 " . ($multiple_samples_per_individual ? "concordance=$con ratio=$ratio" : "ratio=$ratio concordance=$con");
        
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
