#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES EGA_UPLOAD_JAR EGA_DROPBOX EGA_DROPBOX_PW)]);
    use TestPipelines;
}

my $output_dir         = get_output_dir('upload_ega_pipeline');
my $ega_jar            = $ENV{EGA_UPLOAD_JAR};
my $ega_dropbox        = $ENV{EGA_DROPBOX};
my $ega_dropbox_passwd = $ENV{EGA_DROPBOX_PW};

my $testdir = file(qw(t data))->absolute->stringify;
system("rm $testdir/upload_log*");

ok my $pipeline = VRPipe::Pipeline->create(name => 'upload_ega'), 'able to get the upload_ega pipeline';
my @s_names;

foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(ega_upload);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my upload_ega pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data hs_chr20.bam.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        cleanup              => 0,
        'ega_upload_jar'     => "$ega_jar",
        'ega_dropbox'        => "$ega_dropbox",
        'ega_dropbox_passwd' => "$ega_dropbox_passwd",
        'allow_dups'         => 1,
    }
);
ok handle_pipeline(), 'pipeline ran ok';

#*** needs proper tests...

done_testing;
