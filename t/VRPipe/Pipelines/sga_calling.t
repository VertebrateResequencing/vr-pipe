#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
                    required_exe => [qw(sga)]);
    use TestPipelines;
}

ok my $preprocess_pipeline = VRPipe::Pipeline->get(name => 'sga_prepare_fastq'), 'able to get the sga_prepare_fastq pipeline';
ok my $sga_pipeline = VRPipe::Pipeline->get(name => 'sga_variant_calling'), 'able to get the sga_variant_calling pipeline';

my $calling_dir = get_output_dir('sga_calling_test');

my $original_ref_fa = VRPipe::File->get(path => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz))->absolute);
my $ref_fa = file($calling_dir, 'human_g1k_v37.chr11.chr20.fa')->absolute->stringify;
my $new_ref_fa = VRPipe::File->get(path => $ref_fa);
my $oh = $original_ref_fa->openr;
my $nh = $new_ref_fa->openw;
while (<$oh>) {
    print $nh $_;
}
close($oh);
close($nh);

VRPipe::PipelineSetup->get(name => 'sga prepare fastq test',
                           datasource => VRPipe::DataSource->get(type => 'fofn_with_metadata',
                                                                 method => 'grouped_by_metadata',
                                                                 source => file(qw(t data sga_calling_datasource.fofn))->absolute->stringify,
                                                                 options => { metadata_keys => 'sample|platform' },),
                           output_root => $calling_dir,
                           pipeline => $preprocess_pipeline,
                           options => { split_bam_only => '^(11|20)$',
                                        split_bam_include_mate => 1,
                                        sga_exe => 'sga-dindel',
                                        cleanup => 0,
                                      });

VRPipe::PipelineSetup->get(name => 'sga calling test',
                           datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                 method => 'group_by_metadata',
                                                                 source => '1[4]',
                                                                 options => { metadata_keys => 'split_sequence' }),
                           output_root => $calling_dir,
                           pipeline => $sga_pipeline,
                           options => { reference_fasta => $ref_fa,
                                        sga_exe => 'sga-dindel',
                                        cleanup => 0 });

ok handle_pipeline(), 'pipeline ran ok';
ok handle_pipeline(), 'pipeline ran ok 2';

finish;