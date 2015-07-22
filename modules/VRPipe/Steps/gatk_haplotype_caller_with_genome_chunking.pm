
=head1 NAME

VRPipe::Steps::gatk_haplotype_caller_with_genome_chunking - a step

=head1 DESCRIPTION

*** more documentation to come

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

# java -jar GenomeAnalysisTK.jar
#      -T HaplotypeCaller
#      -R reference/human_g1k_v37.fasta
#      -I sample1.bam [-I sample2.bam ...] \
#      --dbsnp dbSNP.vcf \
#      -stand_call_conf [50.0] \
#      -stand_emit_conf 10.0 \
#      [-L targets.interval_list]
#      -o output.raw.snps.indels.vcf

class VRPipe::Steps::gatk_haplotype_caller_with_genome_chunking extends VRPipe::Steps::gatk_haplotype_caller with VRPipe::StepGenomeChunkingRole {
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{bam_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta, caller => 'GATK_HaplotypeCaller' };
            my $options = $self->handle_override_options($vcf_meta);
            $self->handle_standard_options($options);
            
            my $reference_fasta = $options->{reference_fasta};
            my $haplotyper_opts = $options->{haplotype_caller_options};
            my $minimum_records = $options->{minimum_records};
            my $tabix           = $options->{tabix_exe};
            my $avx             = $options->{avx_lsf_requirement_string};
            
            if ($haplotyper_opts =~ /$reference_fasta|-I |--input_file|-o | --output|HaplotypeCaller/) {
                $self->throw("haplotype_caller_options should not include the reference, input or output options or HaplotypeCaller task command");
            }
            
            my $file_list_id;
            if (@{ $self->inputs->{bam_files} } > 1) {
                $file_list_id = VRPipe::FileList->create(files => $self->inputs->{bam_files})->id;
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => qq[java \$jvm_args -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R \$reference_fasta -I \$bams_list -o \$vcf_file $haplotyper_opts -L \$region]
                )
            );
            
            my ($cpus) = $haplotyper_opts =~ m/-nct\s*(\d+)/;
            unless ($cpus) {
                ($cpus) = $haplotyper_opts =~ m/--num_cpu_threads_per_data_thread\s*(\d+)/;
            }
            my $req = $self->new_requirements(memory => 8000, time => 1, $cpus ? (cpus => $cpus) : (), $avx ? (custom => { lsf => $avx }) : ());
            my $basename = 'gatk_haplotype.vcf.gz';
            
            my $chunks = $self->chunks();
            foreach my $chunk (@$chunks) {
                next if (exists $$vcf_meta{chrom} && $$vcf_meta{chrom} ne $chunk->{chrom});
                my $chrom          = $chunk->{chrom};
                my $from           = $chunk->{from};
                my $to             = $chunk->{to};
                my $chunk_opts     = "$haplotyper_opts -L ${chrom}:${from}-${to}";
                my $chunk_basename = "${chrom}_${from}-${to}.$basename";
                my $chunk_meta     = { %$vcf_meta, %$chunk };
                
                my $bams_list_path;
                if ($file_list_id) {
                    $bams_list_path = $self->output_file(basename => "$chunk_basename.bams.list", type => 'txt', temporary => 1)->path;
                }
                else {
                    $bams_list_path = $self->inputs->{bam_files}->[0]->path;
                }
                
                my $vcf_file       = $self->output_file(output_key => 'gatk_hc_vcf_files',       basename => $chunk_basename,       type => 'vcf', metadata => $chunk_meta);
                my $vcf_index_file = $self->output_file(output_key => 'gatk_hc_vcf_index_files', basename => "$chunk_basename.tbi", type => 'bin', metadata => $chunk_meta);
                my $vcf_path       = $vcf_file->path;
                
                $chunk_opts .= qq[ && $tabix -f -p vcf $vcf_path] if ($chunk_opts =~ m/--disable_auto_index_creation_and_locking_when_reading_rods/);
                my $cmd = 'q[' . $self->gatk_prefix($req->memory) . qq[ -T HaplotypeCaller -R $reference_fasta -I $bams_list_path -o $vcf_path $chunk_opts] . ']';
                $cmd .= qq[, input_file_list => $file_list_id] if $file_list_id;
                my $this_cmd = "use VRPipe::Steps::gatk_haplotype_caller; VRPipe::Steps::gatk_haplotype_caller->genotype_and_check($cmd, minimum_records => $minimum_records);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file, $vcf_index_file] });
            }
        };
    }
    
    method outputs_definition {
        return {
            gatk_hc_vcf_files       => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a single vcf file'),
            gatk_hc_vcf_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'a single vcf file'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run GATK HaplotypeCaller for one or more BAMs, generating one compressed VCF per set of BAMs. Sites list can be provided";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
