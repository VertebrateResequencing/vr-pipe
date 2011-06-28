#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Path::Class qw(file dir);

BEGIN {
    use Test::Most tests => 9;
    
    use_ok('VRPipe::Persistent');
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}


# we'll set up 2 simple pipelines: the first has a single step but multiple
# data elements; the second has multiple elements but a single data element.
# all of the following set up stuff should have been tested in PersistentReal.t
# so we don't bother with tests unless we actually get to the manager
my $scheduler = VRPipe::Scheduler->get();
my $output_root = dir($scheduler->output_root, 'manager_test_output');
$scheduler->remove_tree($output_root);
$scheduler->make_path($output_root);

# pipeline 1
my $pipeline1_output_def = VRPipe::FileDefinition->get(name => 'element_output_file',
                                                       type => 'txt',
                                                       output_sub => sub  { my ($self, $step) = @_;
                                                                            return $step->data_element->result.'.o'; });

my $single_step = VRPipe::Step->get(name => 'element_outputter',
                                    inputs_definition => { },
                                    body_sub => sub { my $self = shift;
                                                      my $element_name = $self->data_element->result;
                                                      my $ofile = $self->outputs->{the_only_output}->path;
                                                      $self->dispatch([qq{sleep 5; echo "output for $element_name" > $ofile}, $self->new_requirements(memory => 60, time => 1)]);
                                                      return 0; },
                                    post_process_sub => sub { return 1 },
                                    outputs_definition => { the_only_output => $pipeline1_output_def },
                                    description => 'outputs the data element result to a file');

my $five_element_datasource = VRPipe::DataSource->get(type => 'list', method => 'all', source => file(qw(t data datasource.fivelist)));

my $single_step_pipeline = VRPipe::Pipeline->get(name => 'single_step_pipeline', description => 'simple test pipeline with only a single step');
$single_step_pipeline->add_step($single_step);

my $first_pipeline_output_dir = dir($output_root, 'single_step_pipeline');
my $first_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $five_element_datasource, output_root => $first_pipeline_output_dir, pipeline => $single_step_pipeline);

my @first_output_files;
foreach my $element_num (1..5) {
    push(@first_output_files, file($first_pipeline_output_dir, "fed_result_$element_num.o"));
}

# pipeline 2
my $second_pipeline_output_dir = dir($output_root, 'multi_step_pipeline');
$scheduler->make_path($second_pipeline_output_dir);
my $input1_path = file($second_pipeline_output_dir, 'input1.txt');
open(my $fh, '>', $input1_path) or die "Could not write to $input1_path\n";
print $fh "input1_line1\ninput2_line2\n";
close($fh);
my $input1_file = VRPipe::File->get(path => $input1_path, type => 'txt');
my $output1_path = file($second_pipeline_output_dir, 'output1.txt');
my $output1_file = VRPipe::File->get(path => $output1_path, type => 'txt');

my $single_element_datasource = VRPipe::DataSource->get(type => 'list', method => 'all', source => file(qw(t data datasource.onelist)));

my @steps;
$steps[0] = VRPipe::Step->get(name => 'step_1',
                              inputs_definition => { static_input => $input1_file },
                              body_sub => sub { my $ofile = shift->outputs->{step1_output};
                                                my $fh = $ofile->openw();
                                                print $fh "step1output\n";
                                                $ofile->close();
                                                return 1; },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step1_output => $output1_file },
                              description => 'the first step');
$steps[1] = VRPipe::Step->get(name => "step_2",
                              inputs_definition => { step1_output => VRPipe::FileDefinition->get(name => 'step2_input', type => 'txt') },
                              body_sub => sub { my $ofile = shift->outputs->{step2_output};
                                                my $fh = $ofile->openw();
                                                print $fh "step2output\n";
                                                $ofile->close();
                                                return 1; },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step2_output => VRPipe::FileDefinition->get(name => 'step2_output', type => 'txt') }); #*** would be good to have a test that shows if type => 'bam', the pipeline throws since a txt file was created
$steps[2] = VRPipe::Step->get(name => "step_3",
                              inputs_definition => { step2_output => VRPipe::FileDefinition->get(name => 'step3_input', type => 'txt') },
                              body_sub => sub { my $self = shift;
                                                my $ofile = $self->outputs->{step3_output}->path;
                                                $self->dispatch(["sleep 5; echo step3output > $ofile", $self->new_requirements(memory => 50, time => 1)]);
                                                $self->dispatch(["sleep 4;", $self->new_requirements(memory => 50, time => 1)]);
                                                $self->dispatch(["sleep 3;", $self->new_requirements(memory => 50, time => 1)]);
                                                return 0; },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step3_output => VRPipe::FileDefinition->get(name => 'step3_output', type => 'txt') });
$steps[3] = VRPipe::Step->get(name => "step_4",
                              inputs_definition => { step3_output => VRPipe::FileDefinition->get(name => 'step4_input', type => 'txt') },
                              body_sub => sub { my $ofile = shift->outputs->{step4_output};
                                                my $fh = $ofile->openw();
                                                print $fh "step4output\n";
                                                $ofile->close();
                                                return 1; },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step4_output => VRPipe::FileDefinition->get(name => 'step4_output', output_sub => sub { 'step4_basename.o' }, type => 'txt') });
$steps[4] = VRPipe::Step->get(name => "step_5",
                              inputs_definition => { step4_output => VRPipe::FileDefinition->get(name => 'step5_input', type => 'txt') },
                              body_sub => sub { my $ofile = shift->outputs->{step5_output};
                                                my $fh = $ofile->openw();
                                                print $fh "step5output\n";
                                                $ofile->close();
                                                return 1; },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step5_output => VRPipe::FileDefinition->get(name => 'step5_output', type => 'txt') });

my $multi_step_pipeline = VRPipe::Pipeline->get(name => 'multi_step_pipeline', description => 'simple test pipeline with five steps');
foreach my $step (@steps) {
    $multi_step_pipeline->add_step($step);
}
my $second_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps2', datasource => $single_element_datasource, output_root => $second_pipeline_output_dir, pipeline => $multi_step_pipeline);

my @second_output_files = ($output1_file->path);
foreach my $step_num (2..5) {
    my $ofile = file($second_pipeline_output_dir, $step_num == 4 ? 'step4_basename.o' : "step${step_num}_output_output.txt");
    push(@second_output_files, $ofile);
}

# setup and use manager
ok my $manager = VRPipe::Manager->get(), 'can get a manager with no args';
is $manager->id, 1, 'manager always has an id of 1';
my @manager_setups = $manager->setups;
is @manager_setups, 2, 'setups() returns the correct number of PipelineSetups';
@manager_setups = $manager->setups(pipeline_name => 'multi_step_pipeline');
is @manager_setups, 1, 'setups() returns the correct number of PipelineSetups when pipeline_name supplied';
$manager->trigger;
my $give_up = 100;
while (! $manager->handle_submissions) {
    last if $give_up-- <= 0;
    sleep(1);
}
$give_up = 200;
while (! $manager->trigger) {
    last if $give_up-- <= 0;
    $manager->handle_submissions;
    sleep(1);
}

# check output files for both pipelines
my $all_created = 1;
foreach my $ofile (@first_output_files, @second_output_files) {
    unless (-s $ofile) {
        warn "$ofile is missing\n";
        $all_created = 0;
    }
}
is $all_created, 1, 'multi-step pipeline completed via Manager';

# now lets create a pipeline using a pre-written step, where we'll test that a
# step can work with both a datasource input and the outputs of a previous step
my $prewritten_step_pipeline = VRPipe::Pipeline->get(name => 'md5_pipeline', description => 'simple test pipeline with two steps that create md5 files');
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));

my $fofn_datasource = VRPipe::DataSource->get(type => 'list', method => 'all', source => file(qw(t data datasource.fofn)));
my $prewritten_step_pipeline_output_dir = dir($output_root, 'md5_pipeline');
my $md5_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps3', datasource => $fofn_datasource, output_root => $prewritten_step_pipeline_output_dir, pipeline => $prewritten_step_pipeline);

my @md5_output_files = (file($prewritten_step_pipeline_output_dir, 'file.bam.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.cat.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.txt.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.bam.md5.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.cat.md5.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.txt.md5.md5'));

$give_up = 200;
while (! $manager->trigger) {
    last if $give_up-- <= 0;
    $manager->handle_submissions;
    sleep(1);
}
$all_created = 1;
foreach my $ofile (@md5_output_files) {
    unless (-s $ofile) {
        warn "$ofile is missing\n";
        $all_created = 0;
    }
}
is $all_created, 1, 'all md5 files were created via Manager';

my @md5s;
foreach my $file (file(qw(t data file.bam))->absolute,
                  file(qw(t data file.cat))->absolute,
                  file(qw(t data file.txt))->absolute,
                  file($prewritten_step_pipeline_output_dir, 'file.bam.md5'),
                  file($prewritten_step_pipeline_output_dir, 'file.cat.md5'),
                  file($prewritten_step_pipeline_output_dir, 'file.txt.md5')) {
    push(@md5s, VRPipe::File->get(path => $file)->md5);
}
is_deeply [@md5s], [qw(21efc0b1cc21390f4dcc97795227cdf4 2f8545684149f81e26af90dec0c6869c eb8fa3ffb310ce9a18617210572168ec bc5e5541094af2ddf06d8d6a3ef6e101 3f07e8796553d8dbebdc55d59febeab6 cb683a2796e39c24338cfc23ded5a299)], 'md5s were all set in db';

#*** want to test a datasource that is the outputs of a step of a given (different) pipelinesetup

done_testing;
exit;