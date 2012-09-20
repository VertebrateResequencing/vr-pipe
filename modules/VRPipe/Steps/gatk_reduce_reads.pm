
=head1 NAME

VRPipe::Steps::gatk_reduce_reads - a step

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

# java -Xmx4g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T ReduceReads \
#   -I myData.bam \
#   -o myData.reduced.bam

class VRPipe::Steps::gatk_reduce_reads extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{ $self->$orig }, reduce_reads_options => VRPipe::StepOption->create(description => 'command line options for GATK ReduceReads', optional => 1, default_value => '-l INFO --disable_bam_indexing'), };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more coordinate-sorted bam files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref         = $options->{reference_fasta};
            my $reduce_opts = $options->{reduce_reads_options};
            if ($reduce_opts =~ /$ref|-I |--input_file|-o | --output|ReduceReads/) {
                $self->throw("reduce_reads_options should not include the reference, input files, output files or ReduceReads task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T ReduceReads -R $reference_fasta -I $bam_file -o $reduced_bam_file ' . $reduce_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $reduced_base = $bam->basename;
                $reduced_base =~ s/bam$/reduced.bam/;
                my $reduced_bam_file = $self->output_file(
                    output_key => 'reduced_bam_files',
                    basename   => $reduced_base,
                    type       => 'bam',
                    metadata   => $bam->metadata
                );
                
                my $temp_dir = $options->{tmp_dir} || $reduced_bam_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T ReduceReads -R $ref -I ] . $bam->path . qq[ -o ] . $reduced_bam_file->path . qq[ $reduce_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_reduce_reads', 'reduce_and_check', [$this_cmd, $req, { output_files => [$reduced_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            reduced_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with recalibrated quality scores; OQ tag holds the original quality scores',
                metadata    => { reads => 'Number of reads in the reduced BAM file' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Applies GATK ReduceReads to one or more BAM files which reduces the BAM file using read based compression that keeps only essential information for variant calling";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method reduce_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /-I (\S+) -o (\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $output_reads = $out_file->num_records;
        
        my $bam_type = VRPipe::FileType->create('bam', { file => $out_path });
        
        if ($bam_type->check_magic) {
            $out_file->add_metadata({ reads => $output_reads });
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because output file $out_path does not have the correct magic");
        }
    }
}

1;
