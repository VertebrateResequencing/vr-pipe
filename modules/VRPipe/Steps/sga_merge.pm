=head1 NAME

VRPipe::Steps::sga_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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


# Usage: sga merge [OPTION] ... READS1 READS2
# Merge the sequence files READS1, READS2 into a single file/index
# 
#   -v, --verbose                        display verbose output
#       --help                           display this help and exit
#   -t, --threads=NUM                    use NUM threads to merge the indices (default: 1)
#   -p, --prefix=PREFIX                  write final index to files starting with PREFIX (the default is to concatenate the input filenames)
#   -r, --remove                         remove the original BWT, SAI and reads files after the merge
#   -g, --gap-array=N                    use N bits of storage for each element of the gap array. Acceptable values are 4,8,16 or 32. Lower
#                                        values can substantially reduce the amount of memory required at the cost of less predictable memory usage.
#                                        When this value is set to 32, the memory requirement is essentially deterministic and requires ~5N bytes where
#                                        N is the size of the FM-index of READS2.
#                                        The default value is 4.
#       --no-sequence                    Suppress merging of the sequence files. Use this option when merging the index(es) separate e.g. in parallel
#       --no-forward                     Suppress merging of the forward index. Use this option when merging the index(es) separate e.g. in parallel
#       --no-reverse                     Suppress merging of the reverse index. Use this option when merging the index(es) separate e.g. in parallel

use VRPipe::Base;

class VRPipe::Steps::sga_merge with VRPipe::StepRole {
    method options_definition {
        return { sga_merge_options => VRPipe::StepOption->get(description => 'options to sga index', optional => 1, default_value => '-a bwtsw'),
                 sga_exe => VRPipe::StepOption->get(description => 'path to your sga executable', optional => 1, default_value => 'sga') };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', min_files => 2, max_files => 2, description => 'fastq files to be indexed'),
                 bwt_files => VRPipe::StepIODefinition->get(type => 'bin', min_files => 2, max_files => 2, description => 'fastq files to be indexed') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            
            my $sga_exe = $options->{sga_exe};
            my $sga_opts = $options->{sga_merge_options};
            if ($sga_opts =~ /$ref|merge/) {
                $self->throw("sga_merge_options should not include the merge subcommand");
            }
            my $cmd = $sga_exe.' merge '.$sga_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($bwa_exe, '^Version: (.+)$'), summary => 'sga merge '.$sga_opts.' $fastq_file_1 $fastq_file_2'));
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $fq (@{$self->inputs->{fastq_files}}) {
                my $prefix = $fq->basename;
                $prefix = s/\.(fq|fasta)(\.gz)?//;
                my $bwt_file = $self->output_file(output_key => 'sga_index_binary_files', basename => "$prefix.bwt", type => 'bin', metadata => $fq->metadata);
                $self->dispatch([$cmd, $req, {output_files => [$bwt_file]}]);
            }
        };
    }
    method outputs_definition {
        return { sga_index_binary_files => VRPipe::StepIODefinition->get(type => 'bin', description => 'the files produced by sga index', max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Build the BWT and FM-index for a set of reads";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
