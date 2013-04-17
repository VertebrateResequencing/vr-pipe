
=head1 NAME

VRPipe::Steps::apply_bam_spatial_filter - a step

=head1 DESCRIPTION

Runs spatial_filter to apply previously generated filter file to a set of bams.
Generates either filtered bams (default), or files with problematic reads
marked with the QC fail bit 0x200.

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

class VRPipe::Steps::apply_bam_spatial_filter extends VRPipe::Steps::bamcheck {
    around options_definition {
        my $options = $self->$orig;
        return {
            %{ $self->$orig },
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
                type        => 'any',
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
            my $bamcheck_exe       = $options->{bamcheck_exe};
            my $bc_opts            = VRPipe::Steps::bamcheck->get_bamcheck_options($options);
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            # Get the filter file lane metadata, there should only be one per lane
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
                my $out_log  = $self->output_file(output_key => 'filter_logs', basename => "$basename.log", type => 'txt');
                
                my $qc_fail = $mark_qcfail ? '-f ' : '';
                
                $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'spatial_filter', version => VRPipe::StepCmdSummary->determine_version("$spatial_filter_exe -v", '^spatial_filter: Version v(\S+)'), summary => 'spatial_filter -a -F $filter_path ' . $qc_fail . '$bam_path '));
                
                my $cmd = qq[use VRPipe::Steps::apply_bam_spatial_filter; VRPipe::Steps::apply_bam_spatial_filter->apply_filter(bam_path => q[$bam_path], out_path => q[$out_path], spatial_filter_exe => '$spatial_filter_exe', qc_fail => '$qc_fail', filter_path => '$filter_path', bamcheck_exe => '$bamcheck_exe', bc_opts => '$bc_opts');];
                
                if ($mark_qcfail) {
                    $self->dispatch_vrpipecode($cmd, $req, { output_files => [$out_file, $out_log] });
                }
                else {
                    my $check_file = $self->output_file(output_key => 'bamcheck_files', basename => $out_file->basename . '.bamcheck', type => 'txt');
                    $self->dispatch_vrpipecode($cmd, $req, { output_files => [$out_file, $out_log, $check_file] });
                }
            }
        
        };
    }
    
    method outputs_definition {
        return {
            filtered_bams => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => 'bam with spatially correlated errors either QC failed or filtered out',
                max_files   => -1
            ),
            filter_logs => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'spatial filter log',
                max_files   => -1
            ),
            bamcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'the output of bamcheck on a bam',
                min_files   => 0,
                max_files   => -1,
                metadata    => {
                    source_bam => 'path to the bam file this bamcheck file was created from',
                    lane       => 'lane name (a unique identifer for this sequencing run, aka read group)'
                }
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
        return 0;          # meaning unlimited
    }
    
    method apply_filter (ClassName|Object $self: Str|File :$bam_path!, Str|File :$out_path!, Str :$spatial_filter_exe!, Str :$qc_fail!, Str|File :$filter_path!, Str :$bamcheck_exe!, Str :$bc_opts! ) {
        my $cmd_line = "$spatial_filter_exe -a -F $filter_path $qc_fail $bam_path > $out_path 2>$out_path.log";
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $in_file  = VRPipe::File->get(path => $bam_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        my $log_file = VRPipe::File->get(path => "$out_path.log");
        
        if ($qc_fail) {    # not null, check record count, no need to rebuild metadata
            
            $out_file->update_stats_from_disc(retries => 3);
            my $expected_recs = $in_file->num_records;
            my $actual_recs   = $out_file->num_records;
            
            unless ($actual_recs == $expected_recs) {
                $self->throw("cmd [$cmd_line] failed, $actual_recs recs in the output bam, $expected_recs recs in the original bam");
            }
        }
        else {
            VRPipe::Steps::bamcheck->stats_from_bamcheck("$bamcheck_exe $bc_opts $out_path > $out_path.bamcheck");
            
            # Reconcile read counts before and after filtering with log
            my $reads_before = $in_file->metadata->{reads};
            
            $out_file->reselect_values_from_db; # force metadata refresh
            my $reads_after = $out_file->metadata->{reads};
            
            $log_file->update_stats_from_disc(retries => 3);
            my $logh = $log_file->openr;
            my ($processed, $removed);
            while (my $line = <$logh>) {
                if ($line =~ /Processed\s+(\d+) traces/) {
                    $processed = $1;
                }
                if ($line =~ /Removed\s+(\d+) traces/) {
                    $removed = $1;
                }
            }
            $self->throw("Cannot parse $out_path.log") unless defined $processed && defined $removed;
            $self->throw("Faied to reconcile read counts: before - after [$reads_before - $reads_after] not equal to filtered $removed") unless $reads_before - $reads_after == $removed;
        
        }
        return 1;
    }

}

1;
