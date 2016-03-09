
=head1 NAME

VRPipe::Steps::freebayes_with_genome_chunking - a step

=head1 DESCRIPTION

...

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::freebayes_with_genome_chunking with VRPipe::StepGenomeChunkingRole {
    method options_definition {
        return {
            reference_fasta        => VRPipe::StepOption->create(description => 'absolute path to reference genome fasta'),
            samtools_exe           => VRPipe::StepOption->create(description => 'path to samtools executable for merge', optional => 1, default_value => 'samtools'),
            samtools_merge_options => VRPipe::StepOption->create(description => 'options for samtools merge excluding the -R and -b options', optional => 1, default_value => '-u'),
            freebayes_exe     => VRPipe::StepOption->create(description => 'path to freebayes executable',                                                                                                                                                                                          default_value => 'freebayes'),
            freebayes_options => VRPipe::StepOption->create(description => 'freebayes options',                                                                                                                                                                                                     optional      => 1, default_value => ''),
            bcftools_exe      => VRPipe::StepOption->create(description => 'path to bcftools executable',                                                                                                                                                                                           optional      => 1, default_value => 'bcftools'),
            post_processing   => VRPipe::StepOption->create(description => 'vcf postprocessing command. Must include either $output_vcf or $output_bcf placeholders; can optionally include $bcftools and $reference_fasta as placeholders that will be filled in by the values for those options', default_value => '$bcftools view -Oz > $output_vcf'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files       => VRPipe::StepIODefinition->create(type => 'aln', max_files => -1, description => '1 or more BAM or CRAM files from which to call variants'),
            bam_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => '1 or more index files for BAM or CRAM files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $bcf_meta = $self->common_metadata($self->inputs->{bam_files});
            $bcf_meta = { %$bcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($bcf_meta);
            
            my $samtools       = $options->{samtools_exe};
            my $merge_opts     = $options->{samtools_merge_options};
            my $freebayes      = $options->{freebayes_exe};
            my $freebayes_opts = $options->{freebayes_options};
            my $bcftools       = $options->{bcftools_exe};
            my $post_proc      = $options->{post_processing};
            
            if ($post_proc !~ /\$output_vcf/ && $post_proc !~ /\$output_bcf/) {
                $self->throw("post_proc must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            
            my $reference_fasta = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $reference_fasta->is_absolute;
            
            my $file_list_id;
            if (@{ $self->inputs->{bam_files} } > 1) {
                $file_list_id = VRPipe::FileList->create(files => $self->inputs->{bam_files})->id;
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools merge $merge_opts -R \$region -b \$bams_list -"
                )
            );
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'freebayes',
                    version => VRPipe::StepCmdSummary->determine_version(qq[$freebayes --version], '^version:\s+(\S+)$'),
                    summary => "freebayes -c -f \$reference_fasta $freebayes_opts"
                )
            );
            if ($post_proc =~ m/\$bcftools/) {
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => VRPipe::StepCmdSummary->determine_version($bcftools, '^Version: (.+)$'),
                    summary => "$post_proc"
                );
            }
            
            $post_proc =~ s/\$reference_fasta/$reference_fasta/g;
            $post_proc =~ s/\$bcftools/$bcftools/g;
            
            # define output file
            my $type   = $post_proc =~ /\$output_vcf/ ? 'vcf'    : 'bcf';
            my $suffix = $type eq 'vcf'               ? 'vcf.gz' : 'bcf';
            my $idx    = $type eq 'vcf'               ? 'tbi'    : 'csi';
            my $basename = "freebayes.$suffix";
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            
            my $chunks = $self->chunks();
            foreach my $chunk (@$chunks) {
                my $chrom          = $chunk->{chrom};
                my $from           = $chunk->{from};
                my $to             = $chunk->{to};
                my $chunk_opts     = "$merge_opts -R ${chrom}:${from}-${to}";
                my $chunk_basename = "${chrom}_${from}-${to}.$basename";
                my $chunk_meta     = { %$bcf_meta, %$chunk };
                
                my $bams_list_path;
                if ($file_list_id) {
                    $bams_list_path = $self->output_file(basename => "$chunk_basename.list", type => 'txt', temporary => 1)->path;
                    $bams_list_path = "-b $bams_list_path -";
                }
                else {
                    $bams_list_path = '- ' . $self->inputs->{bam_files}->[0]->path;
                }
                my $output_file  = $self->output_file(output_key => 'freebayes_variant_files',       basename => $chunk_basename,        type => $type, metadata => $chunk_meta);
                my $output_index = $self->output_file(output_key => 'freebayes_variant_index_files', basename => "$chunk_basename.$idx", type => $idx,  metadata => $chunk_meta);
                my $output_path  = $output_file->path;
                my $this_post_proc = $post_proc;
                $this_post_proc =~ s/\$output_[vb]cf/$output_path/g;
                
                my $cmd = qq[$samtools merge $chunk_opts $bams_list_path | $freebayes -c -f $reference_fasta $freebayes_opts | $this_post_proc];
                my $idx_opt = $type eq 'vcf' ? '-ft' : '-f';
                $cmd = "q[($cmd) && $bcftools index $idx_opt $output_path]";
                $cmd .= qq[, input_file_list => $file_list_id] if $file_list_id;
                my $this_cmd = "use VRPipe::Steps::freebayes_with_genome_chunking; VRPipe::Steps::freebayes_with_genome_chunking->freebayes_and_check($cmd, output_path => q[$output_path]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$output_file, $output_index] });
            }
        };
    }
    
    method outputs_definition {
        return {
            freebayes_variant_files       => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => 'a VCF or BCF file per-chunk called across all input BAM/CRAM files'),
            freebayes_variant_index_files => VRPipe::StepIODefinition->create(type => 'idx', max_files => -1, description => 'a .bcf.csi file for each set of input bam files')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run samtools merge and freebayes on one or more BAM or CRAM files. Split calling into multiple jobs across the genome generating one VCF/BCF file per chunk.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method freebayes_and_check (ClassName|Object $self: Str $cmd_line!, File|Str :$output_path!, Int :$input_file_list?) {
        my $output_file = VRPipe::File->get(path => $output_path);
        
        if ($input_file_list) {
            my @input_files = VRPipe::FileList->get(id => $input_file_list)->files;
            my $fofn = VRPipe::File->get(path => qq[$output_path.list]);
            $fofn->create_fofn(\@input_files);
        }
        
        $output_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        return 1;
    }
}

1;
