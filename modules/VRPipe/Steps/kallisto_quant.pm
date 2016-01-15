
=head1 NAME

VRPipe::Steps::kallisto_quant - a step

=head1 DESCRIPTION

RNA-seq quantification with kallisto.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2016 Genome Research Limited.

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

class VRPipe::Steps::kallisto_quant with VRPipe::StepRole {
    method options_definition {
        return {
            kallisto_exe => VRPipe::StepOption->create(
                description   => 'path to your kallisto executable',
                optional      => 1,
                default_value => 'kallisto'
            ),
            kallisto_quant_options => VRPipe::StepOption->create(
                description   => 'options to kallisto quant, excluding -i, -o and --plaintext',
                optional      => 1,
                default_value => '-b 100 -t 4'
            ),
            kallisto_h5dump => VRPipe::StepOption->create(
                description   => 'set to 1 to convert HDF5-formatted results to plaintext',
                optional      => 1,
                default_value => 0
            ),
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools exe',
                optional      => 1,
                default_value => 'samtools'
            ),
        };
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'single-end or paired-end fastq files',
            ),
            transcriptome_index => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => 1,
                description => 'transcriptome index file generated by kallisto index',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $kallisto_exe  = $options->{kallisto_exe};
            my $kallisto_opts = $options->{kallisto_quant_options};
            if ($kallisto_opts =~ /\s-[iop]/) {
                $self->throw("kallisto_quant_options should not include -i, -o, --plaintext or --pseudobam");
            }
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => "$kallisto_exe",
                    version => VRPipe::StepCmdSummary->determine_version($kallisto_exe, '^version (.+)$'),
                    summary => "kallisto quant -i \$fasta_index -o out_dir $kallisto_opts \$reads_1.fq.gz \$reads_2.fq.gz"
                )
            );
            
            my $transcripts_fastq = $self->inputs->{transcriptome_index}->[0]->path;
            my @fastq_files       = map { $_->path } @{ $self->inputs->{fastq_files} };
            my $meta              = $self->common_metadata($self->inputs->{fastq_files});
            
            my $h5_file   = $self->output_file(output_key => 'h5_file',   basename => "abundance.h5",  type => 'bin', metadata => $meta);
            my $tsv_file  = $self->output_file(output_key => 'tsv_file',  basename => "abundance.tsv", type => 'txt', metadata => $meta);
            my $json_file = $self->output_file(output_key => 'json_file', basename => "run_info.json", type => 'txt', metadata => $meta);
            my @out_files = ($h5_file, $tsv_file, $json_file);
            
            my $this_cmd = "$kallisto_exe quant -i $transcripts_fastq -o " . $h5_file->dir . " $kallisto_opts @fastq_files";
            my ($flag, $cpus) = $kallisto_opts =~ /-(t |threads=)(\S+)/;
            if ($kallisto_opts =~ /--pseudobam/) {
                undef $cpus; #pseudobam is not compatible with running on many threads
                $this_cmd =~ s/-t (\S+)//;
                $this_cmd =~ s/--threads (\S+)//;
                my $output_bam = $self->output_file(output_key => 'bam_file', basename => "pseudo_align.bam", type => 'bam', metadata => $meta);
                $this_cmd .= " | " . $options->{samtools_exe} . " view -Sb - > " . $output_bam->path;
                push(@out_files, $output_bam);
            }
            if ($options->{kallisto_h5dump}) {
                my ($flag, $num_bootstraps) = $kallisto_opts =~ /-(b |bootstrap-samples=)(\S+)/;
                foreach my $b (0 .. ($num_bootstraps - 1)) {
                    my $bs_file = $self->output_file(sub_dir => "h5dump", output_key => 'h5dump_files', basename => "bs_abundance_$b.tsv", type => 'txt');
                    push(@out_files, $bs_file);
                }
                $this_cmd .= " && $kallisto_exe h5dump -o ./ " . $h5_file->path;
            }
            my $req = $self->new_requirements(memory => 5000, time => 1, $cpus ? (cpus => $cpus) : ());
            $self->dispatch_wrapped_cmd('VRPipe::Steps::kallisto_quant', 'quantify', [$this_cmd, $req, { output_files => \@out_files }]);
        };
    }
    
    method outputs_definition {
        return {
            h5_file => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => 1,
                description => 'HDF5 binary file containing bootstrap estimates',
            ),
            tsv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => 1,
                description => 'plain-text file of abundance estimates',
            ),
            json_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => 1,
                description => 'json file containing run information',
            ),
            bam_file => VRPipe::StepIODefinition->create(
                type        => 'bam',
                min_files   => 0,
                max_files   => 1,
                description => 'bam file containing pseudo-alignments'
            ),
            h5dump_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 0,
                max_files   => -1,
                description => 'individual bootstrap results'
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs kallisto to quantify abundances of transcripts in single-end or paired-end reads.";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method quantify (ClassName|Object $self: Str $cmd_line) {
        my ($flag, $trans_idx) = $cmd_line =~ /-(i |index=)(\S+)/;
        my $trans_idx_file = VRPipe::File->get(path => $trans_idx);
        $trans_idx_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
    }

}

1;
