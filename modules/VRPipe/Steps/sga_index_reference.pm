
=head1 NAME

VRPipe::Steps::sga_index_reference - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

# Usage: sga index [OPTION] ... READSFILE
# Index the reads in READSFILE using a suffixarray/bwt
#
#   -v, --verbose                        display verbose output
#       --help                           display this help and exit
#   -a, --algorithm=STR                  BWT construction algorithm. STR must be SAIS (induced copying, the default) or BCR (Bauer-Cox-Rosone)
#                                        SAIS is the default method and works well for all types of input. BCR is a specialized to handle
#                                        large volumes of short (<150bp) reads. If you have a large collection of 100bp reads, use BCR as it
#                                        will be much faster and use less memory.
#   -d, --disk=NUM                       use disk-based BWT construction algorithm. The suffix array/BWT will be constructed
#                                        for batchs of NUM reads at a time. To construct the suffix array of 200 megabases of sequence
#                                        requires ~2GB of memory, set this parameter accordingly.
#   -t, --threads=NUM                    use NUM threads to construct the index (default: 1)
#   -c, --check                          validate that the suffix array/bwt is correct
#   -p, --prefix=PREFIX                  write index to file using PREFIX instead of prefix of READSFILE
#       --no-reverse                     suppress construction of the reverse BWT. Use this option when building the index
#                                        for reads that will be error corrected using the k-mer corrector, which only needs the forward index
#       --no-forward                     suppress construction of the forward BWT. Use this option when building the forward and reverse index separately
#   -g, --gap-array=N                    use N bits of storage for each element of the gap array. Acceptable values are 4,8,16 or 32. Lower
#                                        values can substantially reduce the amount of memory required at the cost of less predictable memory usage.
#                                        When this value is set to 32, the memory requirement is essentially deterministic and requires ~5N bytes where
#                                        N is the size of the FM-index of READS2.
#                                        The default value is 8.

use VRPipe::Base;

class VRPipe::Steps::sga_index_reference with VRPipe::StepRole {
    method options_definition {
        return {
            sga_index_reference_options => VRPipe::StepOption->create(description => 'options to sga index to index the reference fasta file', optional => 1, default_value => '-d 10 -t 4'),
            sga_exe                     => VRPipe::StepOption->create(description => 'path to your sga executable',                            optional => 1, default_value => 'sga')
        };
    }
    
    method inputs_definition {
        return { reference_fasta => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'reference fasta file') },;
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $ref_file = $self->inputs->{reference_fasta}[0];
            
            my $sga_exe  = $options->{sga_exe};
            my $sga_opts = $options->{sga_index_reference_options};
            if ($sga_opts =~ /index/) {
                $self->throw("sga_index_reference_options should not include the reference or index subcommand");
            }
            my $cmd = $sga_exe . ' index ' . $sga_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga index ' . $sga_opts . ' $reference_fasta'));
            $cmd .= ' ' . $ref_file->path;
            
            my ($cpus) = $sga_opts =~ m/-t\s*(\d+)/;
            unless ($cpus) {
                ($cpus) = $sga_opts =~ m/--threads (\d+)/;
            }
            my $prefix = $ref_file->basename;
            $prefix =~ s/\.(fa|fasta)(\.gz)?//;
            unless ($sga_opts =~ '--no-forward') {
                my $bwt_file = $self->output_file(output_key => 'sga_index_binary_files', output_dir => $ref_file->dir->stringify, basename => "$prefix.bwt", type => 'bin');
                my $sai_file = $self->output_file(output_key => 'sga_index_text_files',   output_dir => $ref_file->dir->stringify, basename => "$prefix.sai", type => 'txt');
            }
            unless ($sga_opts =~ '--no-reverse') {
                my $bwt_file = $self->output_file(output_key => 'sga_index_binary_files', output_dir => $ref_file->dir->stringify, basename => "$prefix.rbwt", type => 'bin');
                my $sai_file = $self->output_file(output_key => 'sga_index_text_files',   output_dir => $ref_file->dir->stringify, basename => "$prefix.rsai", type => 'txt');
            }
            $self->dispatch([$cmd, $self->new_requirements(memory => 32000, time => 1, $cpus ? (cpus => $cpus) : ()), { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return {
            sga_index_binary_files => VRPipe::StepIODefinition->create(type => 'bin', description => 'the binary files produced by sga index', max_files => 2),
            sga_index_text_files   => VRPipe::StepIODefinition->create(type => 'txt', description => 'the text files produced by sga index',   max_files => 2)
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Indexes a reference genome fasta file, making it suitable for use in subsequent sga variant calling";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
