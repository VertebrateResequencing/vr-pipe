
=head1 NAME

VRPipe::Steps::mpileup_bcf_with_genome_chunking - a step

=head1 DESCRIPTION

Runs samtools mpileup for one or more BAMs, generating one BCF file per set of
BAMs

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::mpileup_bcf_with_genome_chunking extends VRPipe::Steps::mpileup_bcf with VRPipe::StepGenomeChunkingRole {
    method inputs_definition {
        return {
            bam_files       => VRPipe::StepIODefinition->create(type => 'aln', max_files => -1, description => '1 or more BAM or CRAM files from which to call variants'),
            bam_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => '1 or more index files for BAM or CRAM files'),
            sites_file      => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0,  max_files   => 1, description => 'Optional sites file for calling only at the given sites'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $bcf_meta = $self->common_metadata($self->inputs->{bam_files});
            $bcf_meta = { %$bcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($bcf_meta);
            
            my $samtools     = $options->{samtools_exe};
            my $bcftools     = $options->{bcftools_exe};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            
            my $reference_fasta = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $reference_fasta->is_absolute;
            
            if ($mpileup_opts =~ /-f|-b|$reference_fasta/) {
                $self->throw("samtools_mpileup_options should not include the reference, the -f or the -b options");
            }
            
            if ($self->inputs->{sites_file}) {
                $self->throw("samtools_mpileup_options cannot contain the -l option if a sites_file is an input to this step") if ($mpileup_opts =~ /-l/);
                my $sites_file = $self->inputs->{sites_file}[0];
                $mpileup_opts .= " -l " . $sites_file->path;
            }
            
            my $file_list_id;
            if (@{ $self->inputs->{bam_files} } > 1) {
                $file_list_id = VRPipe::FileList->create(files => $self->inputs->{bam_files})->id;
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools mpileup $mpileup_opts -r \$region -f \$reference_fasta -b \$bams_list > \$bcf_file && bcftools index \$bcf_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $basename = 'mpileup.bcf';
            
            my $chunks = $self->chunks();
            foreach my $chunk (@$chunks) {
                next if (exists $$bcf_meta{chrom} && $$bcf_meta{chrom} ne $chunk->{chrom});
                my $chrom          = $chunk->{chrom};
                my $from           = $chunk->{from};
                my $to             = $chunk->{to};
                my $chunk_opts     = "$mpileup_opts -r ${chrom}:${from}-${to}";
                my $chunk_basename = "${chrom}_${from}-${to}.$basename";
                my $chunk_meta     = { %$bcf_meta, %$chunk };
                
                my $bams_list_path;
                if ($file_list_id) {
                    $bams_list_path = $self->output_file(basename => "$chunk_basename.bams.list", type => 'txt', temporary => 1)->path;
                    $bams_list_path = "-b $bams_list_path";
                }
                else {
                    $bams_list_path = $self->inputs->{bam_files}->[0]->path;
                }
                my $bcf_file  = $self->output_file(output_key => 'mpileup_bcf_files',       basename => $chunk_basename,       type => 'bcf', metadata => $chunk_meta);
                my $bcf_index = $self->output_file(output_key => 'mpileup_bcf_index_files', basename => "$chunk_basename.csi", type => 'csi', metadata => $chunk_meta);
                my $bcf_path  = $bcf_file->path;
                
                my $cmd = qq[q[$samtools mpileup $chunk_opts -f $reference_fasta $bams_list_path > $bcf_path && $bcftools index -f $bcf_path]];
                $cmd .= qq[, input_file_list => $file_list_id] if $file_list_id;
                my $this_cmd = "use VRPipe::Steps::mpileup_bcf; VRPipe::Steps::mpileup_bcf->mpileup_bcf_and_check($cmd);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$bcf_file, $bcf_index] });
            }
        };
    }
    
    method outputs_definition {
        return {
            mpileup_bcf_files       => VRPipe::StepIODefinition->create(type => 'bcf', max_files => -1, description => 'a .bcf file for each set of input bam files'),
            mpileup_bcf_index_files => VRPipe::StepIODefinition->create(type => 'csi', max_files => -1, description => 'a .bcf.csi file for each set of input bam files')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run mpileup for one or more BAM or CRAM files. Split calling into multiple jobs across the genome generating one BCF file per chunk.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
