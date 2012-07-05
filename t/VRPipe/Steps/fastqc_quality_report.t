#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],#require env FASTQC ?
                    required_exe => [qw(fastqc)]);
    use TestPipelines;
    use_ok('VRPipe::Steps::fastqc_quality_report');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('fastqc_quality_report', 'fastq_files');

is_deeply [$step->id, $step->description], [1, 'Produces quality report using fastqc'], 'fastqc_quality_report step created and has correct description';

VRPipe::File->get(path => file(qw(t data 2822_6_1.fastq))->absolute);

my $setup = VRPipe::PipelineSetup->get(name => 'fastqc_quality_report',
                                       datasource => VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data fastqc_report_datasource.fofn))->absolute),
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => { fastqc_exe => 'fastqc' }
#                                       options => { fastqc_exe => '/nfs/users/nfs_n/nw11/bin/fastqc' }  
);

my @output_subdirs = output_subdirs(1);
my @zip;
push(@zip, file(@output_subdirs, '1_fastqc_quality_report', "2822_6_1_fastqc.zip"));
ok handle_pipeline(@zip ), 'fastqc_quality_report pipeline ran ok, generating the expected zip file';
finish;
