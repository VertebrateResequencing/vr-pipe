#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Cwd;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('verify_bamid_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'verify_bamid'), 'able to get the verify_bamid pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(verify_bamid);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

# Generate fofn with metadata in output_dir with absolute bam/vcf paths
my $ofofnwm =  VRPipe::File->create(path => file($output_dir, 'verifybam.fofnwm')->absolute);
my $ofh = $ofofnwm->openw or die;
my $ifofnwm = VRPipe::File->create(path => file(qw(t data verifybam.fofnwm))->absolute);
my $ifh = $ifofnwm->openr or die;
my $ddir = cwd() . '/t/data';
while (<$ifh>) {
    if (/^path/) {
        $ofh->write($_);
        next;
    }
    chomp;
    my ($bam, $vcf) = split(/\t/, $_);
    $bam = file($bam)->basename;
    $vcf = file($vcf)->basename;
    $ofh->write("$ddir/$bam\t$ddir/$vcf\n");
}

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my verify_bamid pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn_with_metadata',
        method => 'all',
        source => "$output_dir/verifybam.fofnwm"
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'verify_bamid_exe'         => '/lustre/scratch106/user/cj5/verifyBamID/verifyBamID/bin/verifyBamID',
        'verify_bamid_opts'         => '--ignoreRG',
        cleanup                    => 0
    }
);

my (@output_files);
my $element_id = 0;
while (<$ifh>) {
    if (/^path/) {
        $ofh->write($_);
        next;
    }
    chomp;
    my ($bam, $vcf) = split(/\t/, $_);
    $bam =~ s/\.bam$//;
    $element_id++;
    my @output_dirs = output_subdirs($element_id);
    push(@output_files, file(@output_dirs, '1_verify_bamid', "${bam}.log"));
    push(@output_files, file(@output_dirs, '1_verify_bamid', "${bam}.selfSM"));
    push(@output_files, file(@output_dirs, '1_verify_bamid', "${bam}.depthSM"));
}

ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
