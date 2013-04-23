
=head1 NAME

VRPipe::Steps::htscmd_genotype_analysis - a step

=head1 DESCRIPTION

The genotype check status is determined by analysing the htscmd gtcheck output file, and computing the ratio of best-match genotype concurrence to next-best concurrence; results are stored as metadata (gtype_analysis) on the bam file.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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
            min_concurrence_ratio => VRPipe::StepOption->create(
                description   => 'Minimum best-to-next-best concurrence ratio to be used when applying genotype check status',
                optional      => 1,
                default_value => 1.05
            ),
            min_sites => VRPipe::StepOption->create(
                description   => 'Minimum number of sites for the data to be used in the genotype analysis',
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
            my $self      = shift;
            my $options   = $self->options;
            my $min_ratio = $options->{min_concurrence_ratio};
            my $min_sites = $options->{min_sites};
            my $req       = $self->new_requirements(memory => 3900, time => 1);
            
            foreach my $gt_file (@{ $self->inputs->{htscmd_gtcheck_files} }) {
                my $source_bam  = $gt_file->metadata->{source_bam};
                my $gt_file_path = $gt_file->path;
                my $cmd = "use VRPipe::Steps::htscmd_genotype_analysis; VRPipe::Steps::htscmd_genotype_analysis->analyse_htscmd_output(gt_file_path => q[$gt_file_path], source_bam => q[$source_bam], min_ratio => q[$min_ratio], min_sites => q[$min_sites]);";
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
    
    method analyse_htscmd_output (ClassName|Object $self: Str|File :$gt_file_path!, Str|File :$source_bam!, Num :$min_ratio!, Num :$min_sites!) {

        my $gt_file     = VRPipe::File->get(path => $gt_file_path);
        my $meta        = $gt_file->metadata;
        my $expected    = $meta->{expected_sample};
        my $bam_file    = VRPipe::File->get(path => $source_bam);
        $bam_file->disconnect;

        my $found_expected = 0;
        my ($gtype1,$score1, $gtype2, $score2);
        
        my $pipe = "grep -v '^#' $gt_file_path| sort -nr -k1 |";  # sort descending on Confidence
        my $fh;
        open($fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");

        while (<$fh>) {
            if ($found_expected && defined($gtype1) && defined($gtype2)) { last; }

            # 0.435266        0.468085        25905095.2      25      NA20544
            my ($concurrence,$uncertainty,$avg_depth,$sites,$sample) = split;

            next if $sites < $min_sites;

            if ($expected && $sample eq $expected) { $found_expected = 1; }

            if (!defined $gtype1) {
                $gtype1 = $sample;
                $score1 = $concurrence;
            }
            elsif (!defined $gtype2) {
                $gtype2 = $sample;
                $score2 = $concurrence;
            }
        }
        close $fh;
        
        if ($expected && !$found_expected) { $expected = 0; }
        my $expected_gtype2 = ($expected eq $gtype2 && $score1 == $score2) ? 1 : 0;

        # ratio = GT1 score / GT2 score, best-match / next-best
        my $ratio = $score2 != 0 ? $score1 / $score2 : $score1 / 1e-6;
        $ratio = sprintf("%0.3f", $ratio);
        my $gt_status;

        if ($expected_gtype2) {
            $gt_status = "status=confirmed expected=$expected found=$gtype2 ratio=$ratio";
        }
        elsif ($ratio < $min_ratio) {
            if ($expected) {
                $gt_status = "status=unconfirmed expected=$expected found=$gtype1 ratio=$ratio";
            }
            else {
                $gt_status = "status=unknown expected=none found=$gtype1 ratio=$ratio";
            }
        }
        elsif (!$expected) {
            $gt_status = "status=candidate expected=none found=$gtype1 ratio=$ratio";
        }
        elsif ($expected eq $gtype1) {
            $gt_status = "status=confirmed expected=$expected found=$gtype1 ratio=$ratio";
        }
        else {
            $gt_status = "status=wrong expected=$expected found=$gtype1 ratio=$ratio";
        }
        
        my $new_meta = { gtype_analysis => $gt_status };
        $bam_file->add_metadata($new_meta, replace_data => 1);
        $gt_file->add_metadata($new_meta, replace_data => 1);
    }
}

1;
