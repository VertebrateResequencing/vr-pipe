#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
}

my $output_dir = get_output_dir('gatk_realign_discovery');
ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_realignment_around_discovered_indels'), 'able to get the bam_realignment_around_discovered_indels pipeline';

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

VRPipe::PipelineSetup->create(name => 'indel realignment',
			   datasource => VRPipe::DataSource->create(type => 'fofn_with_metadata',
			  					 method => 'all',
			 					 options => { },
								 source => file(qw(t data pombe_bam.fofnwm))->absolute->stringify),
			   output_root => $output_dir,
			   pipeline => $pipeline,
			   options => { cleanup => 0, 
			  	        reference_fasta => $ref_fa });

ok handle_pipeline(), 'pipeline ran';

#*** needs proper tests

$output_dir = get_output_dir('gatk_snp_calling_and_filter');
ok $pipeline = VRPipe::Pipeline->create(name => 'gatk_variant_calling_and_filter_vcf'), 'able to get the gatk_variant_calling_and_filter_vcf pipeline';

@s_names = ();
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

@expected_step_names = qw(bam_index gatk_genotype vcf_index gatk_variant_filter vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $ds = VRPipe::DataSource->create(type => 'vrpipe',
				 method => 'all',
				 options => { },
				 source => 'indel realignment[4]');

VRPipe::PipelineSetup->create(name => 'snp call and filter',
			   datasource => $ds,
			   output_root => $output_dir,
			   pipeline => $pipeline,
			   options => { delete_unfiltered_vcfs => 0, 
			  	        reference_fasta => $ref_fa,
					genotyper_opts => '-rf BadCigar -glm SNP --output_mode EMIT_ALL_SITES -stand_call_conf 30.0 -stand_emit_conf 10.0 --annotation DepthOfCoverage --annotation QualByDepth --annotation BaseQualityRankSumTest --annotation MappingQualityRankSumTest --annotation FisherStrand --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation RMSMappingQuality',
					var_filter_opts => '--filterExpression "QUAL <= 150" --filterName "QUALfilt" --filterExpression "ADalt <= 18" --filterName "ADaltfil" --filterExpression "Dels > 0.00" --filterName "Delsfilt" --filterExpression "PLalt >= 41" --filterName "PLaltfilt" --filterExpression "ADpropalt <= 0.37" --filterName "ADpropaltfilt" --filterExpression "QD <= 0.45" --filterName "QDfilt" --filterExpression "MQ <= 5.90" --filterName "MQfilt" --filterExpression "SB > -0.01" --filterName "SBfilt" --missingValuesInExpressionsShouldEvaluateAsFailing' });

$output_dir = get_output_dir('gatk_indel_calling_and_filter');
VRPipe::PipelineSetup->create(name => 'indel call and filter',
			   datasource => $ds,
			   output_root => $output_dir,
			   pipeline => $pipeline,
			   options => { delete_unfiltered_vcfs => 0, 
			  	        reference_fasta => $ref_fa,
					genotyper_opts => '-rf BadCigar -glm INDEL --output_mode EMIT_ALL_SITES -stand_call_conf 30.0 -stand_emit_conf 10.0 --annotation DepthOfCoverage --annotation QualByDepth --annotation BaseQualityRankSumTest --annotation MappingQualityRankSumTest --annotation FisherStrand --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation RMSMappingQuality',
					var_filter_opts => '--filterExpression "QUAL < 894" --filterName "QUALfilt" --filterExpression "ADpropalt <  0.04" --filterName "ADpropaltfilt" --filterExpression "GQ < 42.1" --filterName "QDfilt" --filterExpression "MQ < 50.5" --filterName "MQfilt" --filterExpression "ADref > 56" --filterName "ADaltfil" --filterExpression "SB > -0.01" --filterName "SBfilt" --missingValuesInExpressionsShouldEvaluateAsFailing' });

ok handle_pipeline(), 'pipeline ran';

#*** needs proper tests

done_testing;
exit;
