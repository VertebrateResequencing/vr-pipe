#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('bam_htscmd_genotype_checking');

ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_htscmd_genotype_checking'), 'able to get the bam_htscmd_genotype_checking pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(bam_index vcf_sites mpileup_vcf vcf_index htscmd_gtcheck htscmd_genotype_analysis);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $vcf_geno_source = file(qw(t data hs_chr20.genotypes.vcf.gz));
my $vcfg_dir = dir($output_dir, 'vcf_geno');
$pipeline->make_path($vcfg_dir);
my $vcf_geno = file($vcfg_dir, 'hs_chr20.genotypes.vcf.gz')->stringify;
copy($vcf_geno_source, $vcf_geno);
copy("$vcf_geno_source.tbi", "$vcf_geno.tbi");

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $ds = VRPipe::DataSource->create(
    type   => 'fofn',
    method => 'all',
    source => file(qw(t data hs_chr20.qc.bam.fofn))->absolute
);

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my bam_htscmd_genotype_checking pipeline setup',
    datasource => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta => $ref_fa,
        genotypes_vcf => $vcf_geno,
        cleanup  => 0
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->add_metadata({ sample => 'NA20544' });
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '3_mpileup_vcf', "mpileup.vcf.gz"));
    push(@final_files,  file(@output_subdirs, '5_htscmd_gtcheck', "mpileup.vcf.gz.gtypex"));
}

ok handle_pipeline(@output_files,@final_files), 'pipeline ran and created all expected output files';

my @found_gtype_analysis;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    my $meta =  VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->metadata;
    push (@found_gtype_analysis, $meta->{gtype_analysis});
    
}
my @expected_gtype_analysis = (
          'status=unconfirmed expected=NA20544 found=NA20588 ratio=1.000',
          'status=confirmed expected=NA20544 found=NA20544 ratio=1.053',
          'status=confirmed expected=NA20544 found=NA20544 ratio=1.105',
          'status=confirmed expected=NA20544 found=NA20544 ratio=1.105',
);

is_deeply \@found_gtype_analysis, \@expected_gtype_analysis, 'pipeline generated the expected genotype analysis metadata';

done_testing;
exit;
