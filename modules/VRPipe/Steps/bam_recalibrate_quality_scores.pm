
=head1 NAME

VRPipe::Steps::bam_recalibrate_quality_scores - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2013 Genome Research Limited.

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

# java -Xmx4g -jar GenomeAnalysisTK.jar \
#  -l INFO \
#  -R /path/to/reference.fasta \
#  -I<realigned.bam>   \
#  -T TableRecalibration \
#   -o<realigned.recalibrated.bam>   \
#   -recalFile my_reads.recal_data.csv

class VRPipe::Steps::bam_recalibrate_quality_scores extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{ $self->$orig }, bam_recalibration_options => VRPipe::StepOption->create(description => 'command line options for GATK TableRecalibration', optional => 1, default_value => '-l INFO --disable_bam_indexing'), };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more coordinate-sorted bam files'),
            bai_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => -1,
                description => 'index files for the input bam files'
            ),
            bam_recalibration_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => '1 or more bam recal files from count covariates step', metadata => { source_bam => 'path to the bam file used to create this recalibration file' })
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $recal_opts = $options->{bam_recalibration_options};
            if ($recal_opts =~ /$ref|recalFile|TableRecalibration/) {
                $self->throw("bam_recalibration_options should not include the reference, recalFile option or TableRecalibration task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T TableRecalibration -R $reference_fasta -recalFile $bam_file.recal_data.csv -I $bam_file -o $recalibrated_bam_file ' . $recal_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            foreach my $recal_file (@{ $self->inputs->{bam_recalibration_files} }) {
                my $bam_path   = $recal_file->metadata->{source_bam};
                my $bam        = VRPipe::File->get(path => $bam_path);
                my $bam_base   = $bam->basename;
                my $bam_meta   = $bam->metadata;
                my $recal_base = $bam_base;
                $recal_base =~ s/bam$/recal.bam/;
                my $recal_bam_file = $self->output_file(
                    output_key => 'recalibrated_bam_files',
                    basename   => $recal_base,
                    type       => 'bam',
                    metadata   => $bam_meta
                );
                
                my $temp_dir = $options->{tmp_dir} || $recal_bam_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T TableRecalibration -R $ref -recalFile ] . $recal_file->path . qq[ -I ] . $bam->path . qq[ -o ] . $recal_bam_file->path . qq[ $recal_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_recalibrate_quality_scores', 'recal_and_check', [$this_cmd, $req, { output_files => [$recal_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            recalibrated_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with recalibrated quality scores; OQ tag holds the original quality scores',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Recalibrate quality scores of each mapped base using GATK";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method recal_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($recal_path, $in_path, $out_path) = $cmd_line =~ /-recalFile (\S+) -I (\S+) -o (\S+)/;
        $recal_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $in_path    || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path   || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $recal_file = VRPipe::File->get(path => $recal_path);
        my $in_file    = VRPipe::File->get(path => $in_path);
        my $out_file   = VRPipe::File->get(path => $out_path);
        
        # parse the recal file header.
        # if all of the reads have been skipped,
        # then recalibration is not necessary
        my $fh       = $recal_file->openr;
        my $no_recal = 0;
        while (<$fh>) {
            $no_recal = 1 if (/^# Fraction Skipped 1 \/ NaN bp$/);
            last unless /^#/;
        }
        $recal_file->close;
        
        $in_file->disconnect;
        if ($no_recal) {
            $in_file->copy($out_file);
        }
        else {
            system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        }
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output bam file, yet there were $expected_reads reads in the original bam file");
        }
    }
}

1;
