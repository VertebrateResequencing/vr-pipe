#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 8;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
                    required_exe => [qw(samtools bcftools zcat grep wc)]);
    use TestPipelines;
}

ok my $calling_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_chunked_mpileup_bcf'), 'able to get the snp_calling_chunked_mpileup_bcf pipeline';

my $calling_dir = get_output_dir('wgs_calling_test');

my $ref_fa           = file(qw(t data human_g1k_v37.chr11.chr20.fa.gz))->absolute->stringify;
my $all_samples_file = file(qw(t data wgs_calling.samples))->absolute->stringify;
my $override_file    = file(qw(t data wgs_calling_override_options))->absolute->stringify;

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data wgs_calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'wgs_calling_datasource.fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
print $ofh scalar <$ifh>;
while (<$ifh>) {
    chomp;
    my ($source_path, @meta) = split(/\t/, $_);
    my $source = file($source_path);
    my $dest = file($calling_dir, $source->basename);
    copy($source, $dest);
    print $ofh join("\t", $dest, @meta);
    print $ofh "\n";
}
$orig_fofn_file->close;
$fofn_file->close;

VRPipe::PipelineSetup->create(name       => 'wgs test mapping',
                              datasource => VRPipe::DataSource->create(type    => 'fofn_with_metadata',
                                                                       method  => 'grouped_by_metadata',
                                                                       source  => file($calling_dir, 'wgs_calling_datasource.fofn'),
                                                                       options => { metadata_keys => 'analysis_group|chrom' },),
                              output_root => $calling_dir,
                              pipeline    => $calling_pipeline,
                              options     => {
                                           genomic_region_file    => $ref_fa . '.fai',
                                           reference_fasta        => $ref_fa,
                                           chrom_list             => '11 20',
                                           chunk_size             => 10000000,
                                           chunk_overlap          => 0,
                                           sample_sex_file        => $all_samples_file,
                                           chunk_override_options => $override_file,
                                           cleanup                => 0, });

my $regions = { 11 => ['1-10000000', '10000001-20000000', '20000001-30000000', '30000001-40000000', '40000001-50000000', '50000001-60000000', '60000001-70000000', '70000001-80000000', '80000001-90000000', '90000001-100000000', '100000001-110000000', '110000001-120000000', '120000001-130000000', '130000001-135006516'],
                20 => ['1-10000000', '10000001-20000000', '20000001-30000000', '30000001-40000000', '40000001-50000000', '50000001-60000000', '60000001-63025520'], };

my @intermediate_files;
my @final_files;
my %chroms = ('11' => 1, '20' => 2);
while (my ($chrom, $element_id) = each %chroms) {
    my @output_subdirs = output_subdirs($element_id, 1);
    push(@intermediate_files, file(@output_subdirs, '1_chunk_genomic_region', 'human_g1k_v37.chr11.chr20.fa.gz.fai.10000000.0.txt'));
    foreach my $region (@{ $$regions{$chrom} }) {
        my ($from, $to) = $region =~ m/(\d+)-(\d+)$/;
        push(@intermediate_files, file(@output_subdirs, '3_mpileup_chunked_bcf', qq[$chrom\_$from-$to.bcf]));
        push(@intermediate_files, file(@output_subdirs, '4_sex_to_ploidy',       qq[$chrom\_$from-$to.samples]));
        push(@intermediate_files, file(@output_subdirs, '5_bcf_to_vcf',          qq[$chrom\_$from-$to.vcf.gz]));
        push(@intermediate_files, file(@output_subdirs, '5_bcf_to_vcf',          qq[$chrom\_$from-$to.vcf.gz.tbi]));
    }
    push(@final_files, file(@output_subdirs, '7_vcf_concat', 'merged.vcf.gz'));
    push(@final_files, file(@output_subdirs, '7_vcf_concat', 'merged.vcf.gz.tbi'));
}

ok handle_pipeline(@intermediate_files, @final_files), 'pipeline ran ok and all calling files were created';

is_deeply [VRPipe::StepState->create(pipelinesetup => 1, stepmember => 3, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->create(pipelinesetup => 1, stepmember => 5, dataelement => 1)->cmd_summary->summary], ['samtools mpileup -DSV -C50 -m2 -F0.0005 -d 10000 -g -r $region -f $reference_fasta $bam_files > $bcf_file', 'bcftools view -p 0.99 -vcgN $bcf_file | bgzip -c > $vcf_file'], 'cmd summaries for the major steps were as expected';

# check override options work
# indel calling is turned off for chromosome 20 by options in override file
# check that no indels called on chrom20 - grepping for INDEL should only pick up the header INFO line
my $chr20_vcf = $final_files[2]->absolute->stringify;
is `zcat $chr20_vcf | grep INDEL | wc -l`, "1\n", 'chunk_override_options option worked as expected';

# check final vcfs have metadata
is_deeply [VRPipe::File->create(path => $final_files[0])->metadata, VRPipe::File->create(path => $final_files[2])->metadata],
  [{   chrom          => '11',
       analysis_group => 'low_coverage' },
    {  chrom          => '20',
       analysis_group => 'low_coverage' }],
  'final merged vcfs have correct metadata';

# run annotation pipeline
my $vep_cache = file(qw(t data vep_cache))->absolute->stringify;

ok my $annotation_pipeline = VRPipe::Pipeline->create(name => 'vcf_vep_annotate'), 'able to get the vcf_vep_annotate pipeline';

VRPipe::PipelineSetup->create(name       => 'wgs test annotation',
                              datasource => VRPipe::DataSource->create(type    => 'vrpipe',
                                                                       method  => 'all',
                                                                       source  => 'wgs test mapping[7:concat_vcf]',
                                                                       options => {}),
                              output_root => $calling_dir,
                              pipeline    => $annotation_pipeline,
                              options     => {
                                           vep_options                => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir $vep_cache",
                                           'vcf2consequences_options' => "--grantham",
                                           cleanup                    => 0 });

ok handle_pipeline(), 'annotation pipeline ran okay';

my @annotation_files;
foreach my $element_id (3, 4) {
    my @output_dirs = output_subdirs($element_id, 2);
    push(@annotation_files, file(@output_dirs, '1_vep_analysis',         "merged.vep.txt"));
    push(@annotation_files, file(@output_dirs, '2_vcf_vep_consequences', "merged.conseq.vcf.gz"));
    push(@annotation_files, file(@output_dirs, '2_vcf_vep_consequences', "merged.conseq.vcf.gz.tbi"));
}

ok handle_pipeline(@annotation_files), 'annotation pipeline created all expected output files';



# Other pipelines to be written:

# call on subsets - EUR, ASN

# recall on merged site list

finish;
