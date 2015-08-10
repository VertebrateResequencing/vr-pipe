#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VEP ENSEMBL VEP_CACHE)],
        required_exe => [qw(vcf2consequences_vep tabix bcftools)]
    );
    use TestPipelines;
}

my $gerp_cache = file(qw(t data gerp_cache))->absolute->stringify;

# testing first pipeline vcf_vep_annotate
my $output_dir = get_output_dir('vcf_vep_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'vcf_vep_annotate'), 'able to get the vcf_vep_annotate pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(vep_analysis vcf_vep_consequences vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my vcf_vep_annotate pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data datasource.vcf_fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'vep_options'              => "--symbol --format vcf --force_overwrite --offline --cache --dir $ENV{VEP_CACHE}",
        'vcf2consequences_options' => "--gerp $gerp_cache",
        'vep_exe'                  => $ENV{VEP},
        cleanup                    => 1
    }
);

my (@output_files);
my $element_id = 0;
foreach my $in ('test1', 'test2') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '2_vcf_vep_consequences', "${in}.conseq.vcf.gz"));
    push(@output_files, file(@output_dirs, '2_vcf_vep_consequences', "${in}.conseq.vcf.gz.tbi"));
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

# testing second pipeline vep_annotate.pm
$output_dir = get_output_dir('vep_annotate_pipeline');

ok $pipeline = VRPipe::Pipeline->create(name => 'vep_annotate'), 'able to get the vep_annotate pipeline';
@s_names = ();
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
@expected_step_names = qw(vep_annotate_with_genome_chunking bcftools_concat);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

$test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my vep_annotate pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data imputation.vcf.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'vep_options'      => "--assembly GRCh37 --everything --allele_number --format vcf --force_overwrite --plugin Blosum6 --offline --cache --dir $ENV{VEP_CACHE}",
        'post_vep_options' => "bcftools annotate -c CHROM,FROM,TO,INFO/GERP -a $gerp_cache/gerp_score.20.bed.gz -h $gerp_cache/gerp_score.bed.gz.hdr -Ob -o \$output_bcf",
        'vep_exe'          => $ENV{VEP},
        ensembl_api_paths  => $ENV{ENSEMBL},
        chrom_list         => '20',
        chunk_size         => 300000,
        reference_fasta    => file(qw(t data hs37d5.fa))->absolute->stringify,
        cleanup            => 1
    }
);

$element_id++;
@output_files = ();
push(@output_files, file(output_subdirs($element_id, 2), '2_bcftools_concat', "merged.vcf.gz"));
push(@output_files, file(output_subdirs($element_id, 2), '2_bcftools_concat', "merged.vcf.gz.tbi"));
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
