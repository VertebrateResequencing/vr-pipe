#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
}

my $output_dir = get_output_dir('gatk_realign_discovery');
ok my $pipeline = VRPipe::Pipeline->get(name => 'bam_realignment_around_discovered_indels'), 'able to get the bam_realignment_around_discovered_indels pipeline';

my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(bam_metadata bam_index gatk_target_interval_creator_discovery bam_realignment_around_discovered_indels);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data pombe_ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'pombe_ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

VRPipe::PipelineSetup->get(name => 'indel realignment',
			   datasource => VRPipe::DataSource->get(type => 'fofn_with_metadata',
			  					 method => 'all',
			 					 options => { },
								 source => file(qw(t data pombe_bam.fofnwm))->absolute->stringify),
			   output_root => $output_dir,
			   pipeline => $pipeline,
			   options => { cleanup => 0, 
			  	        reference_fasta => $ref_fa });

ok handle_pipeline(), 'pipeline ran';

#*** needs proper tests, and some variant calling steps will come here as well

done_testing;
exit;
