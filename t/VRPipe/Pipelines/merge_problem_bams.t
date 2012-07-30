#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES PICARD GATK)],
                    required_exe => [qw(samtools)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::picard');
}

my $picard           = VRPipe::Steps::picard->new(picard_path => $ENV{PICARD});
my $picard_version   = $picard->determine_picard_version();
my $samtools_version = VRPipe::StepCmdSummary->determine_version('samtools', '^Version: (.+)$');

ok my $merge_lanes_pipeline = VRPipe::Pipeline->create(name => 'bam_merge_lanes_and_fix_rgs'), 'able to get the bam_merge_lanes_and_fix_rgs pipeline';
ok my $merge_libraries_pipeline = VRPipe::Pipeline->create(name => 'bam_merge'), 'able to get the bam_merge pipeline';
my $output_dir = get_output_dir('merge_lanes_and_fix_rgs');

# The bams in this test are outputs of the improvement pipeline given sanger
# sequenced bams. One pair is unaltered, with mismatching RG in header and
# records, one pair have the RGs removed from the records, and one pair has
# matching header and record RGs.
my $poor_dir = dir(qw(t data poor_rg_bams))->absolute;

# we'll test the readgroup_sm_from_metadata_key option by setting individual
# metadata for one of the bam pairs
foreach my $basename ('7662_6#5.bam', '7662_7#5.bam') {
    my $bam = VRPipe::File->create(path => file($poor_dir, $basename));
    $bam->add_metadata({ individual => 'overriden_sample' });
}

my $ref_fa_source = file($poor_dir, 'chr1_truncated.fa');
my $ref_dir = dir($output_dir, 'ref');
$merge_lanes_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'chr1_truncated.fa')->stringify;
copy($ref_fa_source, $ref_fa);

# since merge_lanes_and_fix_rgs requires a vrpipe datasource currently (since
# only that has grouping on metadata_keys functionality), we start with a simple
# import pipeline
VRPipe::PipelineSetup->create(name       => 'test bam import',
                              datasource => VRPipe::DataSource->create(type    => 'fofn',
                                                                       method  => 'all',
                                                                       source  => file($poor_dir, 'fofn')->stringify,
                                                                       options => {}),
                              output_root => $output_dir,
                              pipeline    => VRPipe::Pipeline->create(name => 'import_bams'),
                              options     => {});

VRPipe::PipelineSetup->create(name       => 'test merge lanes',
                              datasource => VRPipe::DataSource->create(type    => 'vrpipe',
                                                                       method  => 'group_by_metadata',
                                                                       source  => 'test bam import[1]',
                                                                       options => { metadata_keys => 'sample|platform|library' }),
                              output_root => $output_dir,
                              pipeline    => $merge_lanes_pipeline,
                              options     => {
                                           bam_tags_to_strip                     => 'OQ XM XG XO',
                                           readgroup_sm_from_metadata_key        => 'individual',
                                           bam_merge_keep_single_paired_separate => 1,
                                           bam_merge_memory                      => 200,
                                           cleanup                               => 1 });

VRPipe::PipelineSetup->create(name       => 'test merge libraries',
                              datasource => VRPipe::DataSource->create(type    => 'vrpipe',
                                                                       method  => 'group_by_metadata',
                                                                       source  => 'test merge lanes[4:markdup_bam_files]',
                                                                       options => { metadata_keys => 'sample|platform' }),
                              output_root => $output_dir,
                              pipeline    => $merge_libraries_pipeline,
                              options     => {
                                           bam_merge_keep_single_paired_separate => 0,
                                           reference_fasta                       => $ref_fa,
                                           reference_public_url                  => 'file:t/data/poor_rg_bams/chr1_truncated.fa',
                                           bam_merge_memory                      => 200,
                                           split_bam_make_unmapped               => 1,
                                           cleanup                               => 1,
                                           cleanup_inputs                        => 1,
                                           remove_merged_bams                    => 1 });

ok handle_pipeline(), 'pipelines ran ok';





#*** needs proper tests; so far only manually confirmed that all 3 resulting
#    sample bams have the expected RGs, and record RG tags match RGs in header





finish;
