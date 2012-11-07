#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(retroseq.pl)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('retroseq_analysis_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'retroseq_analysis'), 'able to get the retroseq_analysis pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(bam_index retroseq_discover retroseq_call);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

#create refTE file, needs to contain absolute file paths
my $refTE_file = file($output_dir, 'chr20.ref_types.tab');
my $alu_file   = file(qw(t data chr20.Alu.bed))->absolute->stringify;
my $l1hs_file  = file(qw(t data chr20.L1_HS.bed))->absolute->stringify;
open(REF, ">", $refTE_file) or die $!;
print REF "Alu\t${alu_file}\n";
print REF "L1HS\t${l1hs_file}\n";
close(REF);

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my retroseq_analysis pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data hs_chr20.bam.fofn))
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'retroseq_exe'          => 'retroseq.pl',
        'refTEs_param'          => "$refTE_file",
        'retroseq_ref'          => file(qw(t data human_g1k_v37.chr20.fa))->absolute->stringify,
        'retroseq_call_options' => '-hets -reads 1 -depth 100',
        cleanup                 => 0,
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '2_retroseq_discover', "${in}.cand.tab"));
    push(@final_files,  file(@output_dirs, '3_retroseq_call',     "${in}.rseq.PE.vcf"));
}
ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
