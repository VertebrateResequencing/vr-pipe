=head1 NAME

VRPipe::Steps::sga_reference_based_calling - a step

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

# Usage: sga graph-diff [OPTION] --base BASE.fa --variant VARIANT.fa
# Find and report strings only present in the graph of VARIANT when compared to BASE
# 
#       --help                           display this help and exit
#       -v, --verbose                    display verbose output
#       -b, --base=FILE                  the baseline reads are in FILE
#       -r, --variant=FILE               the variant reads are in FILE
#       -o, --outfile=FILE               write the strings found to FILE
#       -k, --kmer=K                     use K as the k-mer size for variant discovery
#       -x, --kmer-threshold=T           only used kmers seen at least T times
#       -y, --max-branches=B             allow the search process to branch B times when 
#                                        searching for the completion of a bubble (default: 0)
#       -t, --threads=NUM                use NUM computation threads

use VRPipe::Base;

class VRPipe::Steps::sga_reference_based_calling with VRPipe::StepRole {
    method options_definition {
        return { sga_graph_diff_options => VRPipe::StepOption->get(description => 'options to sga index to index the reference fasta file', optional => 1, default_value => '--debruijn --low-coverage -k 51'),
                 sga_exe => VRPipe::StepOption->get(description => 'path to your sga executable', optional => 1, default_value => 'sga') };
    }
    method inputs_definition {
        return { sga_indexed_variant_reads => VRPipe::StepIODefinition->get(type => 'fq', max_files => -1, description => 'file with variant reads'),
                 reference_fasta => VRPipe::StepIODefinition->get(type => 'txt', max_files => 1, description => 'reference fasta file') },
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref_file = $self->inputs->{reference_fasta}[0];
            
            my $sga_exe = $options->{sga_exe};
            my $sga_opts = $options->{sga_graph_diff_options};
            if ($sga_opts =~ /graph-diff/) {
                $self->throw("sga_graph_diff_options should not include the reference or graph-diff subcommand");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga graph-diff '.$sga_opts.' --variant $variant_fastq --reference $reference_fasta'));
            
            my $req = $self->new_requirements(memory => 16000, time => 1);
            foreach my $fq (@{$self->inputs->{sga_indexed_variant_reads}}) {
                my $prefix = $fq->basename;
                $prefix =~ s/\.(fq|fastq)(\.gz)?$//;
                my $base_vcf = $self->output_file(output_key => 'sga_base_vcf_files',  basename => qq[$prefix.base.vcf], type => 'vcf', metadata => $fq->metadata);
                my $variant_vcf = $self->output_file(output_key => 'sga_variant_vcf_files',  basename => qq[$prefix.variant.vcf], type => 'vcf', metadata => $fq->metadata);
                my $strings_fasta = $self->output_file(output_key => 'sga_strings_fasta_files',  basename => qq[$prefix.strings.fa], type => 'txt', metadata => $fq->metadata);
                my $cmd = qq[$sga_exe graph-diff $sga_opts -p $prefix --variant ].$fq->path.qq[ --reference ].$ref_file->path;
                $self->dispatch([$cmd, $req, {output_files => [$base_vcf,$variant_vcf,$strings_fasta]}]);
                # $self->dispatch([$cmd, $req, {output_files => [$base_vcf,$variant_vcf]}]);
            }
        };
    }
    method outputs_definition {
        return { sga_base_vcf_files => VRPipe::StepIODefinition->get(type => 'vcf', description => 'variant calls made by sga graph-diff', max_files => -1),
                 sga_variant_vcf_files => VRPipe::StepIODefinition->get(type => 'vcf', description => 'variant calls made by sga graph-diff', max_files => -1),
                 sga_strings_fasta_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'variant calls made by sga graph-diff', max_files => -1, check_existence => 0)
                };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Find and report strings only present in the graph of VARIANT when compared to REFERENCE";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
