
=head1 NAME

VRPipe::Steps::gatk_print_reads_with_bqsr - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

# java -Xmx2g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T PrintReads \
#   -o output.bam \
#   -I input.bam \
#   --BQSR input.recal_data.grp

class VRPipe::Steps::gatk_print_reads_with_bqsr extends VRPipe::Steps::gatk_print_reads {
    around inputs_definition {
        return { %{ $self->$orig }, bam_recalibration_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => '1 or more bam recal files from gatk_base_quality_score_recalibration step', metadata => { source_bam => 'path to the bam used to create this recalibration file' })
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = $options->{reference_fasta};
            
            my $recal_opts = $options->{print_reads_options };
            if ($recal_opts =~ /$ref|-I |--input_file|-o |--output|BQSR|PrintReads/) {
                $self->throw("print_reads_options  should not include the reference, input files, output files, BQSR option or PrintReads task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T PrintReads -R $reference_fasta --BQSR $bam_file.recal_data.grp -I $bam_file -o $recalibrated_bam_file ' . $recal_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            foreach my $recal_file (@{ $self->inputs->{bam_recalibration_files} }) {
                my $bam_path   = $recal_file->metadata->{source_bam};
                my $bam        = VRPipe::File->get(path => $bam_path);
                my $recal_base = $bam->basename;
                $recal_base =~ s/bam$/recal.bam/;
                my $recal_bam_file = $self->output_file(
                    output_key => 'recalibrated_bam_files',
                    basename   => $recal_base,
                    type       => 'bam',
                    metadata   => $bam->metadata
                );
                
                my $temp_dir = $options->{tmp_dir} || $recal_bam_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T PrintReads -R $ref --BQSR ] . $recal_file->path . qq[ -I ] . $bam->path . qq[ -o ] . $recal_bam_file->path . qq[ $recal_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_print_reads_with_bqsr', 'apply_bqsr_and_check', [$this_cmd, $req, { output_files => [$recal_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            recalibrated_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with recalibrated quality scores; OQ tag holds the original quality scores',
                metadata    => { reads => 'Number of reads in the recalibrated BAM file' }
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
    
    method apply_bqsr_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /-I (\S+) -o (\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            $out_file->add_metadata({reads => $actual_reads});
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output bam file, yet there were $expected_reads reads in the original bam file");
        }
    }
}

1;
