
=head1 NAME

VRPipe::Steps::sga_permute_reference - a step

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

# Usage: sga preprocess [OPTION] --quality-scale=STR READS1 READS2 ...
# Prepare READS1, READS2, ... data files for assembly
# If pe-mode is turned on (pe-mode=1) then if a read is discarded its pair will be discarded as well.
#
#       --help                           display this help and exit
#       -v, --verbose                    display verbose output
#       -o, --out=FILE                   write the reads to FILE (default: stdout)
#           --phred64                    the input reads are phred64 scaled. They will be converted to phred33.
#       -p, --pe-mode=INT                0 - do not treat reads as paired (default)
#                                        1 - reads are paired with the first read in READS1 and the second
#                                        read in READS2. The paired reads will be interleaved in the output file
#                                        2 - reads are paired and the records are interleaved within a single file.
#       -q, --quality-trim=INT           perform Heng Li's BWA quality trim algorithm.
#                                        Reads are trimmed according to the formula:
#                                        argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT
#                                        where l is the original read length.
#       -f, --quality-filter=INT         discard the read if it contains more than INT low-quality bases.
#                                        Bases with phred score <= 3 are considered low quality. Default: no filtering.
#                                        The filtering is applied after trimming so bases removed are not counted.
#       -r, --remove-adapter-fwd=STRING
#       -c, --remove-adapter-rev=STRING  Remove the adapter STRING from input reads.
#                                        Do not use this option if you are planning to use the BCR algorithm for indexing.
#       -m, --min-length=INT             discard sequences that are shorter than INT
#                                        this is most useful when used in conjunction with --quality-trim. Default: 40
#       -h, --hard-clip=INT              clip all reads to be length INT. In most cases it is better to use
#                                        the soft clip (quality-trim) option.
#       --permute-ambiguous              Randomly change ambiguous base calls to one of possible bases.
#                                        For example M will be changed to A or C. If this option is not specified, the
#                                        entire read will be discarded.
#       -s, --sample=FLOAT               Randomly sample reads or pairs with acceptance probability FLOAT.
#       --dust                           Perform dust-style filtering of low complexity reads. If you are performing
#                                        de novo genome assembly, you probably do not want this.
#       --dust-threshold=FLOAT           filter out reads that have a dust score higher than FLOAT (default: 4.0).
#                                        This option implies --dust
#       --suffix=SUFFIX                  append SUFFIX to each read ID

use VRPipe::Base;

class VRPipe::Steps::sga_permute_reference with VRPipe::StepRole {
    method options_definition {
        return { sga_permute_reference_options => VRPipe::StepOption->create(description => 'options to sga preprocess to permute ambiguous base calls in a reference file', optional => 1, default_value => '--permute-ambiguous'),
                 sga_exe                       => VRPipe::StepOption->create(description => 'path to your sga executable',                                                   optional => 1, default_value => 'sga'),
                 reference_fasta               => VRPipe::StepOption->create(description => 'Absolute path to reference fasta file') };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $ref     = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $sga_exe  = $options->{sga_exe};
            my $sga_opts = $options->{sga_permute_reference_options};
            if ($sga_opts =~ /$ref|preprocess/) {
                $self->throw("sga_permute_reference_options should not include the reference or preprocess subcommand");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga preprocess ' . $sga_opts . ' $reference_fasta > $permuted_reference_fasta'));
            
            my $basename = $ref->basename;
            $basename =~ s/(fa|fasta)(\.gz)?/permute.fa/;
            my $outfile = $self->output_file(output_key => 'permuted_reference_fasta', output_dir => $ref->dir->stringify, basename => $basename, type => 'txt');
            
            my $cmd = qq[$sga_exe preprocess $sga_opts $ref > ] . $outfile->path;
            $self->dispatch([$cmd, $self->new_requirements(memory => 16000, time => 1), { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return { permuted_reference_fasta => VRPipe::StepIODefinition->create(type => 'txt', description => 'reference fasta file with permuted ambiguous base calls', max_files => 1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Permute ambiguous base calls in a fasta reference file for use in sga variant calling";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
