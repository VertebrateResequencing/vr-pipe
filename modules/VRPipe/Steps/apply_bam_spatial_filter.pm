
=head1 NAME

VRPipe::Steps::apply_bam_spatial_filter - a step

=head1 DESCRIPTION

Runs spatial_filter to apply previously generated filter file to a set of bams.
Generates either filtered bams, or files with problematic reads marked with the
QC fail bit 0x200 (default).

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

class VRPipe::Steps::apply_bam_spatial_filter with VRPipe::StepRole {
    method options_definition {
        return {
            spatial_filter_exe => VRPipe::StepOption->create(description => 'path to spatial_filter executable'),
            mark_qcfail        => VRPipe::StepOption->create(description => 'boolean; option to mark problematic reads as QCFAIL rather than filter out', optional => 1, default_value => 0),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => '1 or more bam files, grouped by lane metadata if creating lane-wide filter from a specific tag (phiX control)',
                max_files   => -1,
            ),
            filter_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'filter file identifying bam regions with spatially correlated errors',
                metadata    => { lane => 'lane from which filter was generated' },
                max_files   => 1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self               = shift;
            my $options            = $self->options;
            my $spatial_filter_exe = $options->{spatial_filter_exe};
            my $mark_qcfail        = $options->{mark_qcfail};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            # Get the filter file lane metadata, there should only be one per bam set
            my $filter_file = $self->inputs->{filter_files}[0];
            my $filter_path = $filter_file->path;
            my $filter_lane = $filter_file->metadata->{lane};
            
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $bam_path = $bam_file->path;
                my $lane     = $bam_file->metadata->{lane};
                unless ($lane eq $filter_lane) {
                    $self->throw("bam lane $lane mismatches filter lane $filter_lane");
                }
                
                my $basename = $bam_file->basename;
                $basename =~ s/\.bam/.spfilt.bam/;
                my $out_file = $self->output_file(output_key => 'filtered_bams', basename => $basename, type => 'bam', metadata => $bam_file->metadata);
                my $out_path = $out_file->path;
                
                my $qc_fail = $mark_qcfail ? '-f ' : '';
                
                my $cmd = "$spatial_filter_exe -a -F $filter_path $qc_fail $bam_path > $out_path";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::apply_bam_spatial_filter', 'apply_filter', [$cmd, $req, { output_files => [$out_file] }]);
            }
        
        };
    }
    
    method outputs_definition {
        return {
            filtered_bams => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => 'bam with spatially correlated errors either Wc fialed ot filtered out',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs spatial_filter to QC fail or filter out bam regions with spatially correlated errors e.g. bubbles";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method apply_filter (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        if ($cmd_line =~ / -f /) {
            my ($in_path, $out_path) = $cmd_line =~ / (\S+) > (\S+)$/;
            
            my $in_file  = VRPipe::File->get(path => $in_path);
            my $out_file = VRPipe::File->get(path => $out_path);
            
            $out_file->update_stats_from_disc(retries => 3);
            my $expected_reads = $in_file->num_records;
            my $actual_reads   = $out_file->num_records;
            
            unless ($actual_reads == $expected_reads) {
                $self->throw("cmd [$cmd_line] failed, $actual_reads reads in the output bam, $expected_reads reads in the original bam");
            }
        }
        return 1;
    }

}

1;
