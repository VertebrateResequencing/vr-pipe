#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES IRODS_TEST_ROOT)],
        required_exe => [qw(imkdir iput imeta ils)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_to_irods_pipeline');

my $irods_root=$ENV{IRODS_TEST_ROOT};
foreach my $vcf ("${irods_root}/7/4/b/6/e6315beeeb4c76cdd2c1405e2a4/test2.vcf.gz","${irods_root}/e/a/4/9/a03cb8c09ab6c1e21c82bae0152/test1.vcf.gz" ) {
    system("irm","$vcf") unless system("ils $vcf >/dev/null 2>&1");
    system("irm","$vcf.tbi") unless system("ils $vcf.tbi >/dev/null 2>&1");
}

ok my $pipeline = VRPipe::Pipeline->create(name => 'vcf_to_irods'), 'able to get the vcf_to_irods pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(vcf_metadata  vcf_to_irods);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my vcf_to_irods pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data datasource.vcf_fofn))
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        'irods_root' => "$irods_root",
        'study' => "TEST_STUDY",
        'release' => "REL-2012-11-09",
        cleanup => 1
    }
);

ok handle_pipeline(), 'pipeline ran ok';

my @sample_meta;
foreach my $vcf ("${irods_root}/7/4/b/6/e6315beeeb4c76cdd2c1405e2a4/test2.vcf.gz","${irods_root}/e/a/4/9/a03cb8c09ab6c1e21c82bae0152/test1.vcf.gz" ) {
    my $sm = `imeta ls -d $vcf sample | grep value`;
    my (undef,$s) = split(/ /,$sm);
    chomp($s);
    push (@sample_meta,$s);
    system("irm","$vcf");
    system("irm","$vcf.tbi");
}
my @samples=("UK10K_SCOOP5013919","UK10K_SCOOP5013918");
is_deeply \@samples, \@sample_meta, 'the vcfs were setup on irods with the correct sample metadata';

done_testing;
exit;
