#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(vcf-filter vcf-isec vcf-annotate vcf-stats variant_effect_predictor.pl vcf2consequences_vep tabix)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_filter_merge_and_vep_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'vcf_filter_merge_and_vep_annotate'), 'able to get the vcf_filter_merge_and_vep_annotate pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(vcf_multi_filter vcf_index vcf_merge vcf_index vcf_annotate vcf_index vep_analysis vcf_vep_consequences vcf_index vcf_stats);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $filter_opt_file_1 = file(qw(t data uk10k_gatk_20110715.filter))->absolute->stringify;
my $filter_opt_file_2 = file(qw(t data uk10k_mpileup_20110715.filter))->absolute->stringify;

my $annot_file        = file(qw(t data g1k_dbsnp132_annot.tab.gz))->absolute->stringify;
my $annot_desc_file   = file(qw(t data g1k_dbsnp132_annot_desc.txt))->absolute->stringify;
my $annot_2_file      = file(qw(t data annots-rsIDs-AFs.2011-10-05.tab.gz))->absolute->stringify;
my $annot_2_desc_file = file(qw(t data annots-rsIDs-AFs.2011-10-05.tab.gz.desc))->absolute->stringify;
my $vep_cache         = file(qw(t data vep_cache))->absolute->stringify;
my $gerp_cache        = file(qw(t data gerp_cache))->absolute->stringify;

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my vcf_filter_merge_and_vep_annotate pipeline setup',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        options => { delimiter => "\t" },
        source  => file(qw(t data datasource.vcfs))
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'vcf-annotate_options'     => "-a $annot_file -d $annot_desc_file -c CHROM,FROM,REF,ALT,-,-,INFO/KGPilot123,INFO/dbSNP ",
        'vcf-annotate_2_options'   => "-a $annot_2_file -d $annot_2_desc_file -c CHROM,POS,ID,REF,ALT,INFO/AF_AFR,INFO/AF_AMR,INFO/AF_ASN,INFO/AF_EUR,INFO/AF_MAX ",
        'vcf-filter_programs'      => "vcf-filter#vcf-filter",
        'vcf-filter_files'         => "$filter_opt_file_1#$filter_opt_file_2",
        'vep_options'              => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir $vep_cache",
        'vcf2consequences_options' => "--gerp $gerp_cache",
        'vcf-stats_options'        => "-f FILTER",
        cleanup                    => 0
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $in ('test1', 'test2') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '3_vcf_merge',            "${in}.filt.${in}a.filt.merged.vcf.gz"));
    push(@output_files, file(@output_subdirs, '3_vcf_merge',            "${in}.filt.${in}a.filt.merged.vcf.gz.tbi"));
    push(@output_files, file(@output_subdirs, '5_vcf_annotate',         "${in}.filt.${in}a.filt.merged.annot.vcf.gz"));
    push(@output_files, file(@output_subdirs, '5_vcf_annotate',         "${in}.filt.${in}a.filt.merged.annot.vcf.gz.tbi"));
    push(@final_files,  file(@output_subdirs, '7_vep_analysis',         "${in}.filt.${in}a.filt.merged.annot.conseq.vep.txt"));
    push(@final_files,  file(@output_subdirs, '8_vcf_vep_consequences', "${in}.filt.${in}a.filt.merged.annot.conseq.vcf.gz"));
    push(@final_files,  file(@output_subdirs, '8_vcf_vep_consequences', "${in}.filt.${in}a.filt.merged.annot.conseq.vcf.gz.tbi"));
    push(@final_files,  file(@output_subdirs, '10_vcf_stats',           "${in}.filt.${in}a.filt.merged.annot.conseq.vcf.gz.stats"));
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
