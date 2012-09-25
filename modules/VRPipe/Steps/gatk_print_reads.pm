
=head1 NAME

VRPipe::Steps::gatk_print_reads - a step

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
#   --read_filter MappingQualityZero

class VRPipe::Steps::gatk_print_reads extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{ $self->$orig }, print_reads_options => VRPipe::StepOption->create(description => 'command line options for GATK PrintReads', optional => 1, default_value => '-l INFO --disable_bam_indexing'), };
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
            
            my $ref = $options->{reference_fasta};
            
            my $print_opts = $options->{print_reads_options};
            if ($print_opts =~ /$ref|-I |--input_file|-o |--output|PrintReads/) {
                $self->throw("print_reads_options should not include the reference, input files, output files or PrintReads task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T PrintReads -R $reference_fasta -I $bam_file -o $output_bam_file ' . $print_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $printed_base = $bam->basename;
                $printed_base =~ s/bam$/print.bam/;
                my $printed_bam_file = $self->output_file(
                    output_key => 'printed_bam_files',
                    basename   => $printed_base,
                    type       => 'bam',
                    metadata   => $bam->metadata
                );
                
                my $temp_dir = $options->{tmp_dir} || $printed_bam_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T PrintReads -R $ref -I ] . $bam->path . qq[ -o ] . $printed_bam_file->path . qq[ $print_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_print_reads', 'print_and_check', [$this_cmd, $req, { output_files => [$printed_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            printed_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with recalibrated quality scores; OQ tag holds the original quality scores',
                metadata    => { reads => 'Number of reads in the output BAM file' }
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
    
    method print_and_check (ClassName|Object $self: Str $cmd_line) {
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
        
        # output file does not necessarily have the same number of reads as
        # the input file since it could be filtered. we check that the magic
        # is correct on the output file instead.
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
