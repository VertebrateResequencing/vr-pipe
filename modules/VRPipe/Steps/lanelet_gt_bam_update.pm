
=head1 NAME

VRPipe::Steps::lanelet_gt_bam_update - a step

=head1 DESCRIPTION

Updates bam metadata following lanelet Genotype consistency check using htscmd, setting library_gt = 1 for confirmed bams.

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

class VRPipe::Steps::lanelet_gt_bam_update with VRPipe::StepRole {

    method options_definition {
        return {};
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
            my $self = shift;
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $gt_file (@{ $self->inputs->{htscmd_gtcheck_files} }) {
                my $source_bam  = $gt_file->metadata->{source_bam};
                my $gt_file_path = $gt_file->path;
                my $cmd = "use VRPipe::Steps::lanelet_gt_bam_update; VRPipe::Steps::lanelet_gt_bam_update->analyse_htscmd_output(gt_file_path => q[$gt_file_path], source_bam => q[$source_bam]);";
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
        return "Updates bam metadata following lanelet Genotype consistency check";
    }

    method max_simultaneous {
        return 0;
    }
    
    method analyse_htscmd_output (ClassName|Object $self: Str|File :$gt_file_path!, Str|File :$source_bam!) {

        my $gt_file     = VRPipe::File->get(path => $gt_file_path);
        my $meta        = $gt_file->metadata;
        my $expected    = $meta->{expected_sample};

        my $found_expected = 0;
        my ($gtype1,$score1, $gtype2, $score2);
        
        my $pipe = "grep -v '^#' $gt_file_path| sort -nr -k1 |";  # sort descending on Confidence
        my $fh;
        open($fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");

        while (<$fh>) {
            if ($found_expected && defined($gtype1) && defined($gtype2)) { last; }

            # 0.435266        0.468085        25905095.2      25      NA20544
            my ($concurrence,$uncertainty,$avg_depth,$sites,$sample) = split;

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
        
        if ($expected eq $gtype1 || ($expected eq $gtype2 && $score1 == $score2)) {
            my @bam_files = split ('#',$source_bam);
            foreach my $bam_path (@bam_files) {
                my $bam_file = VRPipe::File->get(path => $bam_path);
                my $original_bam_path = $bam_file->metadata->{original_bam};
                my $original_bam = VRPipe::File->get(path => $original_bam_path);

                # skip update if this is a null bam
                my $meta = $original_bam->metadata;
                next if ($meta->{targeted_bases_mapped} && $meta->{targeted_bases_mapped} == 0) 
                    or $meta->{bases_mapped} == 0;

                $original_bam->add_metadata( {library_gt => 1}, replace_data => 1 );
            }
        }
    }

}

1;
