
=head1 NAME

VRPipe::Steps::bam_processing - a step

=head1 DESCRIPTION

This is a generic bam processing/filtering step. It takes a bam file as input, 
processes it using command(s) given in the options and outputs a bam file.

Example cmd_lines: 1) samtools view -F 0xB00 -b -o $output_bam $input_bam 2)
samtools sort -Obam -T samtools_nsort_tmp $input_bam > $output_bam

Multiple pipes, semi-colons or redirections are allowed as long as one output 
bam is produced for each input bam. It is up to the user to provide a working 
command_line that produces correct outputs. The only requirement is that
strings  $input_bam and $output_bam are both specified in the command line.

Note the step should be used with caution as it does not carry out extra checks
 on the outputs as long as they are in bam format.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::bam_processing with VRPipe::StepRole {
    method options_definition {
        return {
            command_line => VRPipe::StepOption->create(
                description => 'Unix command(s) that will be used to produce the $output_bam from $input_bam (E.g. samtools view -hF 0xB00 $input_bam | samtools sort -Obam -T samtools_nsort_tmp - > $output_bam)',
            ),
            input_metadata_to_keep => VRPipe::StepOption->create(
                description => 'comma-separated list of metadata to be carried over to the processed bams',
                optional    => 1
            ),
            sample_metadata_key => VRPipe::StepOption->create(
                description   => 'metadata key for sample name',
                optional      => 1,
                default_value => 'sample'
            ),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'input bams to be processed',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $cmd_line      = $options->{command_line};
            my $sample_key    = $options->{sample_metadata_key};
            my @metadata_keys = split(',', $options->{input_metadata_to_keep});
            
            if ($cmd_line !~ /\$input_bam/ || $cmd_line !~ /\$output_bam/) {
                $self->throw("cmd_line must contain the strings \$input_bam and \$output_bam");
            }
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => '', version => '', summary => q[$cmd_line]));
            
            my $sub_dir = 'a';
            foreach my $in_bam (@{ $self->inputs->{bam_files} }) {
                my $out_meta = {};
                my $in_meta  = $in_bam->metadata;
                foreach my $key (@metadata_keys) {
                    $out_meta->{$key} = $in_meta->{$key};
                }
                $out_meta->{source_bam} = $in_bam->path;
                my $basename = ($in_meta->{$sample_key}) ? $in_meta->{$sample_key} : "processed";
                my $out_bam = $self->output_file(sub_dir => $sub_dir, output_key => 'bam_files', basename => $basename . ".bam", type => 'bam', metadata => $out_meta);
                ++$sub_dir;
                my $this_cmd    = $cmd_line;
                my $input_path  = $in_bam->path;
                my $output_path = $out_bam->path;
                my $output_dir  = $out_bam->dir;
                $this_cmd =~ s/\$input_bam/$input_path/;
                $this_cmd =~ s/\$output_bam/$output_path/;
                my $cmd_line = "cd $output_dir; " . $this_cmd;
                my $req = $self->new_requirements(memory => 500, time => 1);
                $self->dispatch([$cmd_line, $req, { output_files => [$out_bam] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'output bam file produced by ad-hoc command',
                metadata    => { source_bam => 'original bam from which this file was generated' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs the input bams through command line(s) given by the user";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
