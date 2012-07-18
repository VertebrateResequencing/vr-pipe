#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
                    required_exe => [qw(samtools bcftools)]);
    use TestPipelines;
}

my $output_dir = get_output_dir('snp_calling_mpileup_bcf_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'snp_calling_mpileup_bcf'), 'able to get the snp_calling_mpileup_bcf pipeline';
my @s_names;

foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(mpileup_bcf bcf_to_vcf);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my snp_calling_mpileup_bcf pipeline setup',
    datasource => VRPipe::DataSource->create(type    => 'delimited',
                                             method  => 'all_columns',
                                             options => { delimiter => "\t" },
                                             source  => file(qw(t data hs_chr20.bam.fofn))),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        cleanup       => 0,
        interval_list => file(qw(t data hs_chr20.invervals.bed))->absolute->stringify,
        #samtools_mpileup_options => '-C50 -aug -r 20:1-70000',
        samtools_mpileup_options => '-C50 -aug',
        reference_fasta          => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify, });

my (@output_files, @final_files);
my @files = ('hs_chr20.a.bam', 'hs_chr20.c.bam');
my $element_id = 0;
foreach (@files) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    my $file           = 'mpileup.bcf';
    push(@output_files, file(@output_subdirs, '1_mpileup_bcf', $file));
    $file =~ s/bcf$/vcf.gz/;
    push(@output_files, file(@output_subdirs, '2_bcf_to_vcf', $file));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
