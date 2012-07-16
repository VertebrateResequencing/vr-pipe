#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES TRIMMOMATIC)] #require env TRIMMOMATIC ?
                    #required_exe => [qw(fastqc)]
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::trimmomatic');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('trimmomatic', 'fastq_files');
is_deeply [$step->id, $step->description], [1, 'Step for the Trimmomatic Read Trimmer'], 'trimmomatic step created and has correct description';

# run trimmomatic on a fastq file.
my $setup = VRPipe::PipelineSetup->create(name => 'trimmomatic',
                                       datasource => VRPipe::DataSource->create(
                                                                  type => 'delimited', 
                                                                  method => 'all_columns', 
                                                                  options => { delimiter => "\t" }, 
                                                                  source => file(qw(t data trimmomatic_datasource_pe.fofn))->absolute),   
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => { paired_end => 1 }
);

my @output_subdirs = output_subdirs(1);
my @outfiles;
my $outputfile1 = file(@output_subdirs, '1_trimmomatic', "2822_6_1.paired.trim.fastq");
my $outputfile2 = file(@output_subdirs, '1_trimmomatic', "2822_6_1.unpaired.trim.fastq");
my $outputfile3 = file(@output_subdirs, '1_trimmomatic', "2822_6_2.paired.trim.fastq");
my $outputfile4 = file(@output_subdirs, '1_trimmomatic', "2822_6_2.unpaired.trim.fastq");


my $logfile = file(@output_subdirs, '1_trimmomatic', "trimmomatic.log");
push(@outfiles, $outputfile1,$outputfile2,$outputfile3,$outputfile4, $logfile);
ok handle_pipeline(@outfiles), 'trimmomatic pipeline ran ok, generating the expected trimmed file';

