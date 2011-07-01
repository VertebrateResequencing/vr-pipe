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
                                                                            return $step->data_element->result->{line}.'.o'; });

my $single_step = VRPipe::Step->get(name => 'element_outputter',
                                    inputs_definition => { },
                                    body_sub => sub { my $self = shift;
                                                      my $element_name = $self->data_element->result->{line};
                                                      my $ofile = $self->outputs->{the_only_output}->[0]->path;
                                                      $self->dispatch([qq{sleep 5; echo "output for $element_name" > $ofile}, $self->new_requirements(memory => 60, time => 1)]); },
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
my $multi_step_pipeline = VRPipe::Pipeline->get(name => 'multi_step_pipeline', description => 'simple test pipeline with five steps');

my @steps;
$steps[0] = VRPipe::Step->get(name => 'step_1',
                              inputs_definition => { static_input => $input1_file },
                              body_sub => sub { my $ofile = shift->outputs->{step1_output}->[0];
                                                my $fh = $ofile->openw();
                                                print $fh "step1output\n";
                                                $ofile->close(); },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step1_output => $output1_file },
                              description => 'the first step');
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, before_step_number => 2, adaptor_hash => { step2_input => 'step1_output' });
$steps[1] = VRPipe::Step->get(name => "step_2",
                              inputs_definition => { step2_input => VRPipe::FileDefinition->get(name => 'txt_file', type => 'txt') },
                              body_sub => sub { my $ofile = shift->outputs->{step2_output}->[0];
                                                my $fh = $ofile->openw();
                                                print $fh "step2output\n";
                                                $ofile->close(); },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step2_output => VRPipe::FileDefinition->get(name => 'step2_output', type => 'txt') }); #*** would be good to have a test that shows if type => 'bam', the pipeline throws since a txt file was created
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, before_step_number => 3, adaptor_hash => { step3_input => 'step2_output' });
$steps[2] = VRPipe::Step->get(name => "step_3",
                              inputs_definition => { step3_input => VRPipe::FileDefinition->get(name => 'txt_file', type => 'txt') },
                              body_sub => sub { my $self = shift;
                                                my $ofile = $self->outputs->{step3_output}->[0]->path;
                                                $self->dispatch(["sleep 5; echo step3output > $ofile", $self->new_requirements(memory => 50, time => 1)]);
                                                $self->dispatch(["sleep 4;", $self->new_requirements(memory => 50, time => 1)]);
                                                $self->dispatch(["sleep 3;", $self->new_requirements(memory => 50, time => 1)]); },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step3_output => VRPipe::FileDefinition->get(name => 'step3_output', type => 'txt') });
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, before_step_number => 4, adaptor_hash => { step4_input => 'step3_output' });
$steps[3] = VRPipe::Step->get(name => "step_4",
                              inputs_definition => { step4_input => VRPipe::FileDefinition->get(name => 'txt_file', type => 'txt') },
                              body_sub => sub { my $ofile = shift->outputs->{step4_output}->[0];
                                                my $fh = $ofile->openw();
                                                print $fh "step4output\n";
                                                $ofile->close(); },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step4_output => VRPipe::FileDefinition->get(name => 'step4_o_file', output_sub => sub { 'step4_basename.o' }, type => 'txt') });
VRPipe::StepAdaptor->get(pipeline => $multi_step_pipeline, before_step_number => 5, adaptor_hash => { step5_input => 'step4_output' });
$steps[4] = VRPipe::Step->get(name => "step_5",
                              inputs_definition => { step5_input => VRPipe::FileDefinition->get(name => 'txt_file', type => 'txt') },
                              body_sub => sub { my $ofile = shift->outputs->{step5_output}->[0];
                                                my $fh = $ofile->openw();
                                                print $fh "step5output\n";
                                                $ofile->close(); },
                              post_process_sub => sub { return 1 },
                              outputs_definition => { step5_output => VRPipe::FileDefinition->get(name => 'step5_output', type => 'txt') });

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

is handle_pipeline(@first_output_files, @second_output_files), 1, 'multi-step pipeline completed via Manager';

# now lets create a pipeline using a pre-written step, where we'll test that a
# step can work with both a datasource input and the outputs of a previous step
my $prewritten_step_pipeline = VRPipe::Pipeline->get(name => 'md5_pipeline', description => 'simple test pipeline with two steps that create md5 files');
VRPipe::StepAdaptor->get(pipeline => $prewritten_step_pipeline, before_step_number => 1, adaptor_hash => { md5_file_input => 'data_element' });
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));
VRPipe::StepAdaptor->get(pipeline => $prewritten_step_pipeline, before_step_number => 2, adaptor_hash => { md5_file_input => 'md5_files' });
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));

my $fofn_datasource = VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data datasource.fofn)));
my $prewritten_step_pipeline_output_dir = dir($output_root, 'md5_pipeline');
my $md5_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps3', datasource => $fofn_datasource, output_root => $prewritten_step_pipeline_output_dir, pipeline => $prewritten_step_pipeline);

my @md5_output_files = (file($prewritten_step_pipeline_output_dir, 'file.bam.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.cat.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.txt.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.bam.md5.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.cat.md5.md5'),
                        file($prewritten_step_pipeline_output_dir, 'file.txt.md5.md5'));

is handle_pipeline(@md5_output_files), 1, 'all md5 files were created via Manager';

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

# try out the pre-written mapping pipeline
my $ref_fa = file(qw(t data S_suis_P17.fa));
my $mapping_output_dir = dir($output_root, 'mapping_pipeline');
$scheduler->make_path($mapping_output_dir);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->get(name => 'ps4',
                                                       datasource => VRPipe::DataSource->get(type => 'delimited',
                                                                                             method => 'single_column',
                                                                                             source => file(qw(t data datasource.fastqs)),
                                                                                             options => {delimiter => "\t",
                                                                                                         group_by => 1,
                                                                                                         column => 2}),
                                                       output_root => $mapping_output_dir,
                                                       pipeline => VRPipe::Pipeline->get(name => 'mapping'));

my @mapping_output_files;

#is handle_pipeline(@mapping_output_files), 1, 'all mapping files were created via Manager';

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