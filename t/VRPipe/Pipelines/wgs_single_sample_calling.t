#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK)],
        required_exe => [qw(samtools bcftools vcf-annotate vcf-filter vcf-isec tabix)]
    );
    use TestPipelines;
}

ok my $gatk_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_gatk_unified_genotyper_and_annotate'), 'able to get the snp_calling_gatk_unified_genotyper_and_annotate pipeline';
ok my $mpileup_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_mpileup'), 'able to get the snp_calling_mpileup pipeline';

ok my $filter_pipeline = VRPipe::Pipeline->create(name => 'vcf_annotate'), 'able to get the vcf_annotate pipeline';
ok my $merge_pipeline  = VRPipe::Pipeline->create(name => 'vcf_merge'),    'able to get the vcf_merge pipeline';

my $calling_dir = get_output_dir('wgs_single_sample_calling_test');

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data S_suis_P17.fa))->absolute);
my $ref_fa = VRPipe::File->create(path => file($calling_dir, 'S_suis_P17.fa')->absolute);
$original_ref_fa->copy($ref_fa);

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data wgs_single_sample_calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'wgs_single_sample_calling_datasource.fofn'));
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

my $gatk_filter_file    = file(qw(t data uk10k_gatk_20110715.filter))->absolute->stringify;
my $mpileup_filter_file = file(qw(t data uk10k_mpileup_20110715.filter))->absolute->stringify;

VRPipe::PipelineSetup->create(
    name       => 'gatk single sample calling',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'all',
        source  => $fofn_file->path,
        options => {}
    ),
    output_root => $calling_dir,
    pipeline    => $gatk_pipeline,
    options     => {
        reference_fasta           => $ref_fa->path,
        unified_genotyper_options => '--genotype_likelihoods_model BOTH -stand_call_conf 0.0 -stand_emit_conf 0.0',
        'vcf-annotate_options'    => "-f $gatk_filter_file",
        'vcf-annotate_exe'        => 'vcf-filter',
        cleanup                   => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'mpileup single sample calling',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'all',
        source  => $fofn_file->path,
        options => {}
    ),
    output_root => $calling_dir,
    pipeline    => $mpileup_pipeline,
    options     => {
        reference_fasta       => $ref_fa->path,
        post_calling_vcftools => "vcf-annotate --fill-ICF | vcf-filter -f $mpileup_filter_file",
        cleanup               => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'merge gatk and mpileup',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '1[5]|2[2]',
        options => { metadata_keys => 'sample' }
    ),
    output_root => $calling_dir,
    pipeline    => $merge_pipeline,
    options     => { metadata_priority => 'caller#GATK_UnifiedGenotyper:samtools_mpileup_bcftools' }
);

ok handle_pipeline(), 'calling and merge pipelines ran ok';

my (@output_files, @final_files);
my $element_id = 0;
foreach my $sample ('SAMPLE01', 'SAMPLE02') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 1);
    push(@output_files, file(@output_subdirs, '4_gatk_unified_genotyper', "gatk.vcf.gz"));
    push(@output_files, file(@output_subdirs, '5_vcf_annotate',           "gatk.annot.vcf.gz"));
    push(@output_files, file(@output_subdirs, '5_vcf_annotate',           "gatk.annot.vcf.gz.tbi"));
    
    @output_subdirs = output_subdirs($element_id, 2);
    push(@output_files, file(@output_subdirs, '2_mpileup_vcf', "mpileup.vcf.gz"));
    push(@output_files, file(@output_subdirs, '2_mpileup_vcf', "mpileup.vcf.gz.tbi"));
}

foreach my $sample ('SAMPLE01', 'SAMPLE02') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 3);
    push(@output_files, file(@output_subdirs, '2_vcf_merge', "gatk.annot.mpileup.merged.vcf.gz"));
    push(@output_files, file(@output_subdirs, '2_vcf_merge', "gatk.annot.mpileup.merged.vcf.gz.tbi"));

}
ok handle_pipeline(@output_files), 'pipelines created all expected output files';

finish;

