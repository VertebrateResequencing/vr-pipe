#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Path::Class qw(file dir);

BEGIN {
    use Test::Most tests => 8;
    
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
my $single_step = VRPipe::Step->get(name => 'element_outputter',
                                    inputs_definition => { },
                                    body_sub => sub { my $self = shift;
                                                      my $element_name = $self->data_element->result->{line};
                                                      my $ofile = $self->output_file(output_key => 'the_only_output', basename => "$element_name.o", type => 'txt')->path;
                                                      $self->dispatch([qq{sleep 5; echo "output for $element_name" > $ofile}, $self->new_requirements(memory => 60, time => 1)]); },
                                    outputs_definition => { the_only_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'element_name.o file') },
                                    post_process_sub => sub { return 1 },
                                    description => 'outputs the data element result to a file');

my $five_element_datasource = VRPipe::DataSource->get(type => 'list', method => 'all', source => file(qw(t data datasource.fivelist)));

my $single_step_pipeline = VRPipe::Pipeline->get(name => 'single_step_pipeline', description => 'simple test pipeline with only a single step');
$single_step_pipeline->add_step($single_step);

my $first_pipeline_output_dir = dir($output_root, 'single_step_pipeline');
my $first_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps1', datasource => $five_element_datasource, output_root => $first_pipeline_output_dir, pipeline => $single_step_pipeline);

my @first_output_files;
foreach my $element_num (1..5) {
    push(@first_output_files, file($first_pipeline_output_dir, VRPipe::DataElement->get(datasource => $five_element_datasource->id, result => {line => "fed_result_$element_num"})->id, 'element_outputter', "fed_result_$element_num.o"));
}

# pipeline 2
my $second_pipeline_output_dir = dir($output_root, 'multi_step_pipeline');
$scheduler->make_path($second_pipeline_output_dir);
my $input1_path = file($second_pipeline_output_dir, 'input1.txt');
open(my $fh, '>', $input1_path) or die "Could not write to $input1_path\n";
print $fh "input1_line1\ninput2_line2\n";
close($fh);
my $input1_file = VRPipe::File->get(path => $input1_path, type => 'txt');

my $single_element_datasource = VRPipe::DataSource->get(type => 'list', method => 'all', source => file(qw(t data datasource.onelist)));
my $multi_step_pipeline = VRPipe::Pipeline->get(name => 'multi_step_pipeline', description => 'simple test pipeline with five steps');

my $output1_path = file($second_pipeline_output_dir, VRPipe::DataElement->get(datasource => $single_element_datasource->id, result => {line => "single_element"})->id, "step_1", 'output1.txt');
my $output1_file = VRPipe::File->get(path => $output1_path, type => 'txt');

my @steps;
$steps[0] = VRPipe::Step->get(name => 'step_1',
                              inputs_definition => { static_input => $input1_file },
                              body_sub => sub { my $ofile = shift->output_file(output_key => 'step1_output', basename => 'output1.txt', type => 'txt');
                                                my $fh = $ofile->openw();
                                                print $fh "step1output\n";
                                                $ofile->close(); },
                              outputs_definition => { step1_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step1_output file') },
                              post_process_sub => sub { return 1 },
                              description => 'the first step');
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, to_step => 2, adaptor_hash => { step2_input => { step1_output => 1 } });
$steps[1] = VRPipe::Step->get(name => "step_2",
                              inputs_definition => { step2_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_2 input') },
                              body_sub => sub {
                                                my $ofile = shift->output_file(output_key => 'step2_output', basename => 'step2_output.txt', type => 'txt');
                                                my $fh = $ofile->openw();
                                                print $fh "step2output\n";
                                                $ofile->close();
                                               },
                              outputs_definition => { step2_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step2_output file') },
                              post_process_sub => sub { return 1 });
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, to_step => 3, adaptor_hash => { step3_input => { step2_output => 2 } });
$steps[2] = VRPipe::Step->get(name => "step_3",
                              inputs_definition => { step3_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_3 input') },
                              body_sub => sub { 
                                                my $self = shift;
                                                my $ofile = $self->output_file(output_key => 'step3_output', basename => 'step3_output.txt', type => 'txt')->path;
                                                $self->dispatch(["sleep 5; echo step3output > $ofile", $self->new_requirements(memory => 50, time => 1)]);
                                                $self->dispatch(["sleep 4;", $self->new_requirements(memory => 50, time => 1)]);
                                                $self->dispatch(["sleep 3;", $self->new_requirements(memory => 50, time => 1)]);
                                               },
                              outputs_definition => { step3_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step3_output file') },
                              post_process_sub => sub { return 1 });
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, to_step => 4, adaptor_hash => { step4_input => { step3_output => 3 } });
$steps[3] = VRPipe::Step->get(name => "step_4",
                              inputs_definition => { step4_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_4 input') },
                              body_sub => sub {
                                                my $ofile = shift->output_file(output_key => 'step4_output', basename => 'step4_basename.o', type => 'txt');
                                                my $fh = $ofile->openw();
                                                print $fh "step4output\n";
                                                $ofile->close();
                                               },
                              outputs_definition => { step4_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step4_output file') },
                              post_process_sub => sub { return 1 });
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, to_step => 5, adaptor_hash => { step5_input => { step4_output => 4 } });
$steps[4] = VRPipe::Step->get(name => "step_5",
                              inputs_definition => { step5_input => VRPipe::StepIODefinition->get(type => 'txt', description => 'step_5 input') },
                              body_sub => sub {
                                                my $ofile = shift->output_file(output_key => 'step5_output', basename => 'step5_output.txt', type => 'txt');
                                                my $fh = $ofile->openw();
                                                print $fh "step5output\n";
                                                $ofile->close();
                                               },
                              outputs_definition => { step5_output => VRPipe::StepIODefinition->get(type => 'txt', description => 'step5_output file') },
                              post_process_sub => sub { return 1 });

foreach my $step (@steps) {
    $multi_step_pipeline->add_step($step);
}
my $second_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps2', datasource => $single_element_datasource, output_root => $second_pipeline_output_dir, pipeline => $multi_step_pipeline);

my @second_output_files = ($output1_file->path);
foreach my $step_num (2..5) {
    my $ofile = file($second_pipeline_output_dir, VRPipe::DataElement->get(datasource => $single_element_datasource->id, result => {line => "single_element"})->id, "step_$step_num", $step_num == 4 ? 'step4_basename.o' : "step${step_num}_output.txt");
    push(@second_output_files, $ofile);
}

# setup and use manager
ok my $manager = VRPipe::Manager->get(), 'can get a manager with no args';
is $manager->id, 1, 'manager always has an id of 1';
my @manager_setups = $manager->setups;
is @manager_setups, 2, 'setups() returns the correct number of PipelineSetups';
@manager_setups = $manager->setups(pipeline_name => 'multi_step_pipeline');
is @manager_setups, 1, 'setups() returns the correct number of PipelineSetups when pipeline_name supplied';

is handle_pipeline(@first_output_files, @second_output_files), 1, 'multi-step pipeline completed via Manager';

# now lets create a pipeline using a pre-written step, where we'll test that a
# step can work with both a datasource input and the outputs of a previous step
my $prewritten_step_pipeline = VRPipe::Pipeline->get(name => 'md5_pipeline', description => 'simple test pipeline with two steps that create md5 files');
VRPipe::StepAdaptor->get(pipeline => $prewritten_step_pipeline, to_step => 1, adaptor_hash => { md5_file_input => { data_element => 0 } });
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));
VRPipe::StepAdaptor->get(pipeline => $prewritten_step_pipeline, to_step => 2, adaptor_hash => { md5_file_input => { md5_files => => 1 } });
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));

my $fofn_datasource = VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data datasource.fofn)));
my $prewritten_step_pipeline_output_dir = dir($output_root, 'md5_pipeline');
my $md5_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps3',
                                                   datasource => $fofn_datasource,
                                                   output_root => $prewritten_step_pipeline_output_dir,
                                                   pipeline => $prewritten_step_pipeline,
                                                   options => {md5_files_in_source_dir => 0});

my @md5_output_files = (file($prewritten_step_pipeline_output_dir, 7, 'md5_file_production', 'file.bam.md5'),
                        file($prewritten_step_pipeline_output_dir, 8, 'md5_file_production', 'file.cat.md5'),
                        file($prewritten_step_pipeline_output_dir, 9, 'md5_file_production', 'file.txt.md5'),
                        file($prewritten_step_pipeline_output_dir, 7, 'md5_file_production', 'file.bam.md5.md5'),
                        file($prewritten_step_pipeline_output_dir, 8, 'md5_file_production', 'file.cat.md5.md5'),
                        file($prewritten_step_pipeline_output_dir, 9, 'md5_file_production', 'file.txt.md5.md5'));

is handle_pipeline(@md5_output_files), 1, 'all md5 files were created via Manager';

my @md5s;
foreach my $file (file(qw(t data file.bam))->absolute,
                  file(qw(t data file.cat))->absolute,
                  file(qw(t data file.txt))->absolute,
                  file($prewritten_step_pipeline_output_dir, 7, 'md5_file_production', 'file.bam.md5'),
                  file($prewritten_step_pipeline_output_dir, 8, 'md5_file_production', 'file.cat.md5'),
                  file($prewritten_step_pipeline_output_dir, 9, 'md5_file_production', 'file.txt.md5')) {
    push(@md5s, VRPipe::File->get(path => $file)->md5);
}
is_deeply [@md5s], [qw(21efc0b1cc21390f4dcc97795227cdf4 2f8545684149f81e26af90dec0c6869c eb8fa3ffb310ce9a18617210572168ec bc5e5541094af2ddf06d8d6a3ef6e101 3f07e8796553d8dbebdc55d59febeab6 cb683a2796e39c24338cfc23ded5a299)], 'md5s were all set in db';

#*** want to test a datasource that is the outputs of a step of a given (different) pipelinesetup

done_testing;
exit;

sub handle_pipeline {
    my $give_up = 200;
    while (! $manager->trigger) {
        last if $give_up-- <= 0;
        $manager->handle_submissions;
        sleep(1);
    }
    my $all_created = 1;
    foreach my $ofile (@_) {
        unless (-s $ofile) {
            warn "$ofile is missing\n";
            $all_created = 0;
        }
    }
    return $all_created;
}