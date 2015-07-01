
=head1 NAME

VRPipe::Steps::gatk_combine_gvcfs_with_genome_chunking - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>

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

# java -Xmx2g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T CombineGVCFs \
#   --variant gvcf1.vcf \
#   --variant gvcf2.vcf \
#   -o mergeGvcf.vcf

class VRPipe::Steps::gatk_combine_gvcfs_with_genome_chunking extends VRPipe::Steps::gatk_combine_gvcfs with VRPipe::StepGenomeChunkingRole {
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{gvcf_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($vcf_meta);
            $self->handle_standard_options($options);
            
            my $reference_fasta          = $options->{reference_fasta};
            my $combine_gvcfs_opts       = $options->{combine_gvcfs_options} ? $options->{combine_gvcfs_options} : "";
            my $maximum_gvcfs_to_combine = $options->{maximum_gvcfs_to_combine};
            
            if ($combine_gvcfs_opts =~ /$reference_fasta|-V |--variant|-o | --output|CombineGVCFs/) {
                $self->throw("combine_gvcfs_options should not include the reference, input or output options or CombineGVCFs task command");
            }
            
            my $basename = 'gatk_combineGvcf.vcf.gz';
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T CombineGVCFs -R $reference_fasta --variant $gvcf1 --variant $gvcf2 [...] -o $out_gvcf -L $region ' . $combine_gvcfs_opts
                )
            );
            
            my @input_files = @{ $self->inputs->{gvcf_files} };
            
            my $req = $self->new_requirements(memory => 6000, time => 1);
            
            my $chunks = $self->chunks();
            my $count  = 0;
            while (@input_files) {
                my @files_list = splice @input_files, 0, $maximum_gvcfs_to_combine;
                my @file_inputs = map { "--variant " . $_->path } @files_list;
                foreach my $chunk (@$chunks) {
                    next if (exists $$vcf_meta{chrom} && $$vcf_meta{chrom} ne $chunk->{chrom});
                    my $chrom          = $chunk->{chrom};
                    my $from           = $chunk->{from};
                    my $to             = $chunk->{to};
                    my $chunk_opts     = "$combine_gvcfs_opts -L ${chrom}:${from}-${to}";
                    my $chunk_basename = "${chrom}_${from}-${to}.$basename";
                    my $chunk_meta     = { %$vcf_meta, %$chunk };
                    my $vcf_file       = $self->output_file(output_key => 'combine_gvcf_files', basename => "batch_$count." . $chunk_basename, type => 'vcf', metadata => $chunk_meta);
                    my $vcf_path       = $vcf_file->path;
                    my $vcf_index      = $self->output_file(output_key => 'combine_gvcf_index_files', basename => "batch_$count." . $chunk_basename . ".tbi", type => 'bin', metadata => $chunk_meta);
                    
                    my $cmd      = $self->gatk_prefix($req->memory) . qq[ -T CombineGVCFs -R $reference_fasta @file_inputs -o $vcf_path $chunk_opts];
                    my $this_cmd = "use VRPipe::Steps::gatk_combine_gvcfs; VRPipe::Steps::gatk_combine_gvcfs->combine_and_check(q[$cmd]);";
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file, $vcf_index] });
                    $count++;
                }
            }
        };
    }
    
    method description {
        return "Run GATK CombineGVCFs on equal sized batches of input gVCF files produced by the Haplotype Caller, generating one joint VCF file for each batch";
    }

}

1;
