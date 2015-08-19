
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
                description => 'Unix command(s) that will be used to produce the $output_bam or $output_cram from $input_bam or $input_cram (E.g. samtools view -hF 0xB00 $input_bam | samtools sort -Obam -T samtools_nsort_tmp - > $output_bam)',
            ),
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to samtools executable if using $samtools placeholder in command_line',
                optional      => 1,
                default_value => 'samtools'
            ),
            input_metadata_to_keep => VRPipe::StepOption->create(
                description => 'comma-separated list of metadata to be carried over to the processed files. By default copies all metadata to the output files; set to "-"" to transfer no metadata',
                optional    => 1
            ),
            index_output => VRPipe::StepOption->create(
                description   => 'boolean; index the output BAM or CRAM file when true',
                optional      => 1,
                default_value => 0
            ),
            check_records_vs_input => VRPipe::StepOption->create(
                description   => 'boolean; check the number of input records equals the number of output records when true',
                optional      => 1,
                default_value => 0
            ),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'input BAM or CRAM files to be processed',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $cmd_line      = $options->{command_line};
            my $samtools      = $options->{samtools_exe};
            my $idx_output    = $options->{index_output};
            my $check_records = $options->{check_records_vs_input};
            my @metadata_keys = split(',', $options->{input_metadata_to_keep});
            
            if ($cmd_line !~ /\$input_bam/ && $cmd_line !~ /\$input_cram/ && $cmd_line !~ /\$input_aln/) {
                $self->throw("cmd_line must contain one of the input strings \$input_bam, \$input_cram or \$input_aln");
            }
            if ($cmd_line !~ /\$output_bam/ && $cmd_line !~ /\$output_cram/) {
                $self->throw("cmd_line must contain one of the output strings \$output_bam or \$output_cram");
            }
            if ($cmd_line =~ /\$samtools/) {
                $self->set_cmd_summary(
                    VRPipe::StepCmdSummary->create(
                        exe     => 'samtools',
                        version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                        summary => '$cmd_line'
                    )
                );
                $cmd_line =~ s/\$samtools/$samtools/g;
            }
            else {
                $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => '', version => '', summary => q[$cmd_line]));
            }
            
            my $sub_dir = 'a';
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $in_bam (@{ $self->inputs->{bam_files} }) {
                my $out_meta = {};
                my $in_meta  = $in_bam->metadata;
                if (@metadata_keys) {
                    foreach my $key (@metadata_keys) {
                        next if $key eq '-';
                        $out_meta->{$key} = $in_meta->{$key};
                    }
                }
                else {
                    $out_meta = $in_meta;
                }
                $out_meta->{source_bam} = $in_bam->path;
                my $basename = $in_bam->basename;
                $basename =~ s/\.(cr|b)am$//;
                my $suffix = $cmd_line =~ /\$output_bam/ ? 'bam' : 'cram';
                my $out_bam     = $self->output_file(sub_dir => $sub_dir, output_key => 'processed_bam_files', basename => "$basename.$suffix", type => $suffix, metadata => $out_meta);
                my @out_files   = ($out_bam);
                my $input_path  = $in_bam->path;
                my $output_path = $out_bam->path;
                my $output_dir  = $out_bam->dir;
                my $this_cmd    = $cmd_line;
                $this_cmd =~ s/\$input_(cr|b)am/$input_path/g;
                $this_cmd =~ s/\$output_(cr|b)am/$output_path/g;
                
                if ($idx_output) {
                    my $suffix_idx = $cmd_line =~ /\$output_bam/ ? 'bai' : 'crai';
                    my $out_idx = $self->output_file(sub_dir => $sub_dir, output_key => 'processed_index_files', basename => "$basename.$suffix.$suffix_idx", type => 'bin', metadata => $out_meta);
                    push @out_files, $out_idx;
                }
                $this_cmd = "cd $output_dir; " . $this_cmd;
                
                my $args = qq[q[$this_cmd], input_path => q[$input_path], output_path => q[$output_path]];
                $args .= qq[, check_records_vs_input => 1]  if $check_records;
                $args .= qq[, samtools_exe => q[$samtools]] if $idx_output;
                
                my $this_cmd_line = "use VRPipe::Steps::bam_processing; VRPipe::Steps::bam_processing->process_and_check($args);";
                $self->dispatch_vrpipecode($this_cmd_line, $req, { output_files => \@out_files });
                
                ++$sub_dir;
            }
        };
    }
    
    method outputs_definition {
        return {
            processed_bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'output BAM or CRAM files produced by ad-hoc command',
                metadata    => { source_bam => 'original BAM or CRAM from which this file was generated' }
            ),
            processed_index_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                min_files   => 0,
                max_files   => -1,
                description => 'output BAI or CRAI index files',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs the input BAM or CRAM through command line(s) given by the user";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method process_and_check (ClassName|Object $self: Str $cmd_line, Str|File :$input_path!, Str|File :$output_path!, Bool :$check_records_vs_input?, Str :$samtools_exe?) {
        my $input_file  = VRPipe::File->get(path => $input_path);
        my $output_file = VRPipe::File->get(path => $output_path);
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        if ($samtools_exe) {
            system(qq[$samtools_exe index $output_path]) && $self->throw("failed to index the output file [$samtools_exe index $output_path]");
        }
        
        if ($check_records_vs_input) {
            $output_file->_filetype->check_records_vs_input($input_file, $cmd_line);
        }
    }

}

1;
