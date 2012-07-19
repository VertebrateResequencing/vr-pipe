#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;
use Archive::Zip;
use Archive::Zip::MemberRead;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)], #require env FASTQC ?
                    required_exe => [qw(fastqc)]);
    use TestPipelines;
    use_ok('VRPipe::Steps::fastqc_quality_report');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('fastqc_quality_report', 'fastq_files');

is_deeply [$step->id, $step->description], [1, 'Produces quality report using fastqc'], 'fastqc_quality_report step created and has correct description';

VRPipe::File->create(path => file(qw(t data 2822_6_1.fastq))->absolute);

my $setup = VRPipe::PipelineSetup->create(name        => 'fastqc_quality_report',
                                          datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data fastqc_report_datasource.fofn))->absolute),
                                          output_root => $output_dir,
                                          pipeline    => $pipeline,
                                          options => { fastqc_exe => 'fastqc' });

my @output_subdirs = output_subdirs(1);
my @zip;
my $outputfile = file(@output_subdirs, '1_fastqc_quality_report', "2822_6_1_fastqc.zip");
push(@zip, $outputfile);
ok handle_pipeline(@zip), 'fastqc_quality_report pipeline ran ok, generating the expected zip file';

# -- check content of fastqc_data.txt
# Get output file
my $zip = Archive::Zip->new($outputfile->stringify);

my $fh = Archive::Zip::MemberRead->new($zip, "2822_6_1_fastqc/fastqc_data.txt");
my $outputfile_content = '';
my $line;
while (defined($line = $fh->getline())) {
    $outputfile_content .= $line;
}
$fh->close;

#Get reference file in t/data
my $testfile = 't/data/2822_6_1_fastqc.zip';
$zip = Archive::Zip->new($testfile);
$fh = Archive::Zip::MemberRead->new($zip, "2822_6_1_fastqc/fastqc_data.txt");
my $testfile_content = '';
while (defined($line = $fh->getline())) {
    $testfile_content .= $line;
}
$fh->close;
is($testfile_content, $outputfile_content, 'The fastqc report file is the same as the reference test files');

finish;
