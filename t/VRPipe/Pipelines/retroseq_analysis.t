#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
	);
    use TestPipelines;
	use lib "/software/vertres/lib/all/Restroseq";
}

my $output_dir = get_output_dir('retroseq_analysis_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'retroseq_analysis'), 'able to get the retroseq_analysis pipeline';
my @s_names;
print STDERR "DESC=",$pipeline->description,"\n";

foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(retroseq_discover retroseq_call);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my retroseq_analysis pipeline setup',
		datasource => VRPipe::DataSource->get(type => 'fofn',
			method => 'all',
			source => file(qw(t data hs_chr20.bam.fofn))),
		output_root => $output_dir,
		pipeline => $pipeline,
		options => { 
            'retroseq_exe' => '/software/vertres/bin-external/retroseq.pl',
            'refTEs_param' => '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources/TEs/Repeatmasker/ref_types.tab',
            'exde_param' => '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources/TEs/Repeatmasker/exclude_regions.fofn',
            'retroseq_ref' => '/lustre/scratch105/projects/g1k/ref/phase1/human_g1k_v37.fasta',
            'retroseq_filter' => '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources/TEs/Repeatmasker/ref_types.tab',
            'retroseq_call_options' => '-hets -reads 5 -depth 400',
		    cleanup => 0,
        });

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b') {
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '1_retroseq_discover', "${in}.cand.tab"));
    push(@final_files, file(@output_dirs, '2_retroseq_call', "${in}.rseq.vcf.PE"));
}
use Data::Dumper;print STDERR Dumper(@output_files);
use Data::Dumper;print STDERR Dumper(@final_files);
ok handle_pipeline(@output_files,@final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
