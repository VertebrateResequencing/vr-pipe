#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $output_dir = get_output_dir('test_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'test_pipeline'), 'able to get the test_pipeline pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(test_step_one test_step_two test_step_three test_step_four)], 'the pipeline has the correct steps';

$pipeline = VRPipe::Pipeline->get(name => 'test_pipeline');
@s_names = ();
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(test_step_one test_step_two test_step_three test_step_four)], 'the pipeline has the correct steps after a second retrieval';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my test pipeline setup',
                                                    datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                          method => 'all',
                                                                                          source => file(qw(t data datasource.fofn2))),
                                                    output_root => $output_dir,
                                                    pipeline => $pipeline,
                                                    options => {all_option => 'foo',
                                                                one_option => 50,
                                                                four_option => 'bar'});

my @output_files;
foreach my $in ('file.txt', 'file2.txt', 'file3.txt') {
    foreach my $suffix ('step_one', 'step_one.step_two', 'step_one.step_two.step_three', 'step_one.step_two.step_three.step_four') {
        push(@output_files, file($output_dir, "$in.$suffix"));
    }
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

my $ofile = VRPipe::File->get(path => file($output_dir, 'file3.txt.step_one.step_two.step_three.step_four'));
my $ometa = $ofile->metadata;
my $o2meta = VRPipe::File->get(path => file($output_dir, 'file2.txt.step_one.step_two.step_three.step_four'))->metadata;
is_deeply [$ometa->{one_meta}, $ometa->{two_meta}, $ometa->{three_meta}, $o2meta->{three_meta}, $ometa->{four_meta}], [50, 'body_decided_two_option', 'no_three_meta', 'StepOption_default_decided_three_option', 'bar'], 'metadata of one of the final output files was as expected';

done_testing;
exit;