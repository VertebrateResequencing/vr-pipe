
=head1 NAME

VRPipe::Steps::sga_index - a step

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

class VRPipe::Steps::sga_index with VRPipe::StepRole {
    method options_definition {
        return {
            sga_index_options => VRPipe::StepOption->create(description => 'options to sga index',        optional => 1, default_value => '--no-reverse'),
            sga_exe           => VRPipe::StepOption->create(description => 'path to your sga executable', optional => 1, default_value => 'sga')
        };
    }
    
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->create(type => 'fq', max_files => -1, description => 'sequence files to be indexed', metadata => { bases => 'number of bases in the fastq file' }) };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $sga_exe  = $options->{sga_exe};
            my $sga_opts = $options->{sga_index_options};
            if ($sga_opts =~ /index/) {
                $self->throw("sga_index_options should not include the index subcommand");
            }
            my $cmd = $sga_exe . ' index ' . $sga_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga index ' . $sga_opts . ' $reads_file'));
            
            my ($cpus) = $sga_opts =~ m/-t\s*(\d+)/;
            unless ($cpus) {
                ($cpus) = $sga_opts =~ m/--threads (\d+)/;
            }
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                # estimate 5 bits/base for memory
                my $bases  = $fq->metadata->{bases};
                my $memory = 5000 + int(5 * $bases / (8 * 1024 * 1024) + 0.5);                                   # Mb
                my $req    = $self->new_requirements(memory => $memory, time => 1, $cpus ? (cpus => $cpus) : ());
                
                my $prefix = $fq->basename;
                $prefix =~ s/\.(fq|fastq)(\.gz)?//;
                my @outfiles;
                unless ($sga_opts =~ '--no-forward') {
                    push @outfiles, $self->output_file(output_key => 'sga_index_binary_files', output_dir => $fq->dir->stringify, basename => "$prefix.bwt", type => 'bin');
                    push @outfiles, $self->output_file(output_key => 'sga_index_text_files',   output_dir => $fq->dir->stringify, basename => "$prefix.sai", type => 'txt');
                }
                unless ($sga_opts =~ '--no-reverse') {
                    push @outfiles, $self->output_file(output_key => 'sga_index_binary_files', output_dir => $fq->dir->stringify, basename => "$prefix.rbwt", type => 'bin');
                    push @outfiles, $self->output_file(output_key => 'sga_index_text_files',   output_dir => $fq->dir->stringify, basename => "$prefix.rsai", type => 'txt');
                }
                my $this_cmd = "$cmd " . $fq->path;
                $self->dispatch([$this_cmd, $req, { output_files => \@outfiles }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            sga_index_binary_files => VRPipe::StepIODefinition->create(type => 'bin', description => 'the binary files produced by sga index', max_files => -1),
            sga_index_text_files   => VRPipe::StepIODefinition->create(type => 'txt', description => 'the text files produced by sga index',   max_files => -1)
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Build the BWT and FM-index for a set of reads";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
