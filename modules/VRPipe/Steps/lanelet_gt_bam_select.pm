
=head1 NAME

VRPipe::Steps::lanelet_gt_bam_select - a step

=head1 DESCRIPTION

Checks a set of BAMs and selects which ones should be included in a lanelet Genotype consistency check.
Either selects all bams if a mimimum confirmed GBM (Gb bases mapped) has not been reached, or ('fast' mode) selects a subset of confirmed bams plus all unchecked bams.
The bams must have bases-mapped metadata, and should be grouped by sample, eg by library metadata.

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

class VRPipe::Steps::lanelet_gt_bam_select with VRPipe::StepRole {

    method options_definition {
        return {
            minimum_lib_gbm => VRPipe::StepOption->create(description => 'minimum total Gbm for each bam set required to run this check', optional => 1, default_value => 500),
            minimum_gbm_confirmed => VRPipe::StepOption->create(description => 'minimum total Gbm for confirmed bams below which we need to check all bams', optional => 1, default_value => 500),
        };
    }
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to be GT checked'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options   = $self->options;
            my $minimum_lib_gbm = $options->{minimum_lib_gbm};
            my $minimum_gbm_confirmed = $options->{minimum_gbm_confirmed};

            my $total_bases_mapped = 0;  # count of total bases mapped for bam group
            my $confirmed_bases_mapped = 0;  # count of mapped bases already confirmed

            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta = $bam->metadata;
                my $bases_mapped = $meta->{targeted_bases_mapped} ? $meta->{targeted_bases_mapped} : $meta->{bases_mapped};
                $self->throw("Could not find bases_mapped or targeted_bases_mapped in bam metadata") unless $bases_mapped;

                $total_bases_mapped += $bases_mapped;
                if ($meta->{library_gt}) {
                    $confirmed_bases_mapped += $bases_mapped;
                }
            }

            my $req = $self->new_requirements(memory => 500, time => 1);

            # If we do not have enough data from this library yet, create a null bam just to get the pipeline to finish
            if ($total_bases_mapped / 1000000 < $minimum_lib_gbm) {
                ##$self->warn("Insufficent Gb mapped bases " .  $total_bases_mapped / 1000000 . " below minumum $minimum_lib_gbm");
                foreach my $bam (@{ $self->inputs->{bam_files} }) {
                    my $bam_in_path = $bam->path;
                    my $basename = $bam->basename;
                    my $bam_out = $self->output_file(output_key => 'bam_files', basename => $basename, type => 'bam');
                    my $bam_out_path = $bam_out->path;
                    my $cmd = "samtools view -H -b $bam_in_path > $bam_out_path;";

                    $self->dispatch_wrapped_cmd('VRPipe::Steps::lanelet_gt_bam_select', 'copy_bam', [$cmd, $req, { output_files => [$bam_out] }]);
                }
                return;
            }

            # Take a copy of the bams we want to GT check

            if ($confirmed_bases_mapped / 1000000 < $minimum_gbm_confirmed) {
                # Confirmed below the minimum - check all bams
                foreach my $bam (@{ $self->inputs->{bam_files} }) {

                    my $bam_in_path = $bam->path;
                    my $basename = $bam->basename;
                    my $bam_out = $self->output_file(output_key => 'bam_files', basename => $basename, type => 'bam');
                    my $bam_out_path = $bam_out->path;
                    my $cmd = "cp $bam_in_path $bam_out_path;";

                    $self->dispatch_wrapped_cmd('VRPipe::Steps::lanelet_gt_bam_select', 'copy_bam', [$cmd, $req, { output_files => [$bam_out] }]);
                }
            }
            else {
                # Confirmed above the minimum - check all unconfirmed bams with a subset of confirmed bams
                my $total_confirmed_bases = 0;  # Count of bp for all confirmed bams that will be included in GT check

                foreach my $bam (@{ $self->inputs->{bam_files} }) {

                    my $bam_in_path = $bam->path;
                    my $basename = $bam->basename;
                    my $bam_out = $self->output_file(output_key => 'bam_files', basename => $basename, type => 'bam');
                    my $bam_out_path = $bam_out->path;
                    my $cmd = "cp $bam_in_path $bam_out_path;";

                    my $meta = $bam->metadata;
                    if ($meta->{library_gt}) {
                        if ($total_confirmed_bases < $minimum_gbm_confirmed) {
                            $self->dispatch_wrapped_cmd('VRPipe::Steps::lanelet_gt_bam_select', 'copy_bam', [$cmd, $req, { output_files => [$bam_out] }]);

                            my $bases_mapped = $meta->{targeted_bases_mapped} ? $meta->{targeted_bases_mapped} : $meta->{bases_mapped};
                            $total_confirmed_bases += $bases_mapped;
                        }
                    }
                    else {
                        $self->dispatch_wrapped_cmd('VRPipe::Steps::lanelet_gt_bam_select', 'copy_bam', [$cmd, $req, { output_files => [$bam_out] }]);
                    }
                }
            }
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }

    method description {
        return "Selects which BAMs should be included in a lanelet Genotype consistency check";
    }

    method max_simultaneous {
        return 0;            # meaning unlimited
    }

    method copy_bam (ClassName|Object $self: Str $cmd_line) {

        system($cmd_line) && $self->throw("failed to run [$cmd_line]");

        my $null_file = $cmd_line =~ /^cp/ ? 0 : 1;

        my ($input_path, $output_path);
        if ($null_file) {
            ($input_path, $output_path) = $cmd_line =~ / (\S+) > (\S+);$/;
        }
        else {
            ($input_path, $output_path) = $cmd_line =~ /^cp (\S+) (\S+);$/;
        }

        my $meta = VRPipe::File->get(path => $input_path)->metadata;
        if ($null_file) {
            if ($meta->{targeted_bases_mapped}) {
                $meta->{targeted_bases_mapped} = 0;
            }
            else {
                $meta->{bases_mapped} = 0;
            }
        }

        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->add_metadata({  %$meta, original_bam => $input_path} );
    }
    
    method outputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'bam file for GT check'
            )
        };
    }
}

1;
