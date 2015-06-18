
=head1 NAME

VRPipe::Steps::bam_to_fixed_bam - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk> and Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bam_to_fixed_bam with VRPipe::StepRole {
    method options_definition {
        return {
            reference_fasta               => VRPipe::StepOption->create(description => 'absolute path to genome reference file used to do the mapping'),
            uncompressed_fixed_bam_output => VRPipe::StepOption->create(
                description   => 'A boolean to choose if the output bams from this step should be uncompressed',
                optional      => 1,
                default_value => 1
            ),
            fixed_bam_seq_from_reference => VRPipe::StepOption->create(
                description => 'A boolean to choose whether sequence info is read from the reference -- set for mappers, such as smalt, which don\'t include this in output bam files',
                optional    => 1
            ),
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            )
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'raw bam files from a mapper',
                metadata    => {
                    reads  => 'total number of reads (sequences)',
                    paired => '0=unpaired reads were mapped; 1=paired reads were mapped'
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $ref     = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $samtools         = $options->{samtools_exe};
            my $samtools_version = VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$');
            my $fixmate_extra    = '';
            if ($samtools_version =~ /^(\d+)\.\d+/) {
                if ($1 > 0) {
                    $fixmate_extra = ' -O bam';
                }
            }
            
            my $u = $options->{uncompressed_fixed_bam_output} ? ' -u' : ' -b';
            my $t = $options->{fixed_bam_seq_from_reference} ? ' -T $reference_fasta' : '';
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => $samtools_version,
                    summary => "samtools view -bu \$bam_file$t | samtools sort -n -o - samtools_nsort_tmp | samtools fixmate$fixmate_extra /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp | samtools fillmd$u - \$reference_fasta > \$fixed_bam_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 2900, time => 1);
            foreach my $input (@{ $self->inputs->{bam_files} }) {
                my $basename    = $input->basename;
                my $prefix_base = '.' . $basename;
                
                my $output = $self->output_file(
                    output_key => 'fixed_bam_files',
                    basename   => $basename,
                    type       => 'bam',
                    metadata   => $input->metadata
                );
                
                my $dir         = $output->dir;
                my $input_path  = $input->path;
                my $output_path = $output->path;
                my $nprefix     = file($dir, $prefix_base . '.samtools_nsort_tmp');
                my $cprefix     = file($dir, $prefix_base . '.samtools_csort_tmp');
                
                my $tr = $options->{fixed_bam_seq_from_reference} ? " -T $ref" : '';
                my $this_cmd = "$samtools view -bu $input_path$tr | $samtools sort -n -o - $nprefix | $samtools fixmate$fixmate_extra /dev/stdin /dev/stdout | $samtools sort -o - $cprefix | $samtools fillmd$u - $ref > $output_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_to_fixed_bam', 'fix_and_check', [$this_cmd, $req, { output_files => [$output] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            fixed_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'uncompressed coordinate-sorted bam file(s)',
                metadata    => {
                    reads  => 'total number of reads (sequences)',
                    paired => '0=unpaired reads were mapped; 1=paired reads were mapped'
                }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Turns a bam file into an uncompressed coordinate-sorted bam file with fixed mates and correct NM tag values";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method fix_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($samtools, $input_path, $output_path) = $cmd_line =~ /^(\S+) view -bu (\S+) .+ (\S+)$/;
        $output_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $input  = VRPipe::File->get(path => $input_path);
        my $output = VRPipe::File->get(path => $output_path);
        
        $output->disconnect;
        system('set -o pipefail; ' . $cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $output->_filetype->check_records_vs_input($input, $cmd_line);
    }
}

1;
