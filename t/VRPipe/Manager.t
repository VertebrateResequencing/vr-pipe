#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Path::Class qw(file dir);
use POSIX qw(getgroups);
use VRPipe::Interface::CmdLine;

BEGIN {
    use Test::Most tests => 19;
    use VRPipeTest;
    use TestPipelines;
}

# we'll set up 2 simple pipelines: the first has a single step but multiple
# data elements; the second has multiple steps but a single data element.
# all of the following set up stuff should have been tested in Persistent.t
# so we don't bother with tests unless we actually get to the manager
my $scheduler   = VRPipe::Scheduler->create();
my $output_root = get_output_dir('manager_test_output');

# pipeline 1
my $single_step = VRPipe::Step->create(
    name              => 'element_outputter',
    inputs_definition => {},
    body_sub          => sub {
        my $self         = shift;
        my $element_name = $self->data_element->result->{line};
        my $ofile        = $self->output_file(output_key => 'the_only_output', basename => "$element_name.o", type => 'txt')->path;
        $self->dispatch([qq{sleep 5; echo "output for $element_name" > $ofile}, $self->new_requirements(memory => 60, time => 1)]);
    },
    outputs_definition => { the_only_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'element_name.o file') },
    post_process_sub   => sub               { return 1 },
    description        => 'outputs the data element result to a file'
);

my $five_element_datasource = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.fivelist))->absolute);

my $single_step_pipeline = VRPipe::Pipeline->create(name => 'single_step_pipeline', description => 'simple test pipeline with only a single step');
$single_step_pipeline->add_step($single_step);

my $first_pipeline_output_dir = dir($output_root, 'single_step_pipeline');
my $first_pipelinesetup = VRPipe::PipelineSetup->create(name => 'ps1', datasource => $five_element_datasource, output_root => $first_pipeline_output_dir, pipeline => $single_step_pipeline);

my @first_output_files;
$five_element_datasource->elements;
foreach my $element_num (1 .. 5) {
    push(@first_output_files, file(output_subdirs(VRPipe::DataElement->get(datasource => $five_element_datasource->id, result => { line => "fed_result_$element_num" })->id), '1_element_outputter', "fed_result_$element_num.o"));
}

# pipeline 2
my $second_pipeline_output_dir = dir($output_root, 'multi_step_pipeline');
$scheduler->make_path($second_pipeline_output_dir);
my $input1_path = file($second_pipeline_output_dir, 'input1.txt');
open(my $fh, '>', $input1_path) or die "Could not write to $input1_path\n";
print $fh "input1_line1\ninput2_line2\n";
close($fh);
my $input1_file = VRPipe::File->create(path => $input1_path, type => 'txt');

my $single_element_datasource = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.onelist))->absolute);
my $multi_step_pipeline = VRPipe::Pipeline->create(name => 'multi_step_pipeline', description => 'simple test pipeline with five steps');

my @steps;
$steps[0] = VRPipe::Step->create(
    name              => 'step_1',
    inputs_definition => { static_input => $input1_file },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step1_output', basename => 'output1.txt', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step1output\n";
        $ofile->close();
    },
    outputs_definition => { step1_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step1_output file') },
    post_process_sub   => sub            { return 1 },
    description        => 'the first step'
);
VRPipe::StepAdaptor->create(pipeline => $multi_step_pipeline, to_step => 2, adaptor_hash => { step2_input => { step1_output => 1 } });
$steps[1] = VRPipe::Step->create(
    name              => "step_2",
    inputs_definition => { step2_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_2 input') },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step2_output', basename => 'step2_output.txt', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step2output\n";
        $ofile->close();
    },
    outputs_definition => { step2_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step2_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::StepAdaptor->create(pipeline => $multi_step_pipeline, to_step => 3, adaptor_hash => { step3_input => { step2_output => 2 } });
$steps[2] = VRPipe::Step->create(
    name              => "step_3",
    inputs_definition => { step3_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_3 input') },
    body_sub          => sub {
        my $self = shift;
        my $ofile = $self->output_file(output_key => 'step3_output', basename => 'step3_output.txt', type => 'txt')->path;
        $self->dispatch(["sleep 5; echo step3output > $ofile", $self->new_requirements(memory => 50, time => 1)]);
        $self->dispatch(["sleep 4;",                           $self->new_requirements(memory => 50, time => 1)]);
        $self->dispatch(["sleep 3;",                           $self->new_requirements(memory => 50, time => 1)]);
    },
    outputs_definition => { step3_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step3_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::StepAdaptor->create(pipeline => $multi_step_pipeline, to_step => 4, adaptor_hash => { step4_input => { step3_output => 3 } });
$steps[3] = VRPipe::Step->create(
    name              => "step_4",
    inputs_definition => { step4_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_4 input') },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step4_output', basename => 'step4_basename.o', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step4output\n";
        $ofile->close();
    },
    outputs_definition => { step4_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step4_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::StepAdaptor->create(pipeline => $multi_step_pipeline, to_step => 5, adaptor_hash => { step5_input => { step4_output => 4 } });
$steps[4] = VRPipe::Step->create(
    name              => "step_5",
    inputs_definition => { step5_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'step_5 input') },
    body_sub          => sub {
        my $ofile = shift->output_file(output_key => 'step5_output', basename => 'step5_output.txt', type => 'txt');
        my $fh = $ofile->openw();
        print $fh "step5output\n";
        $ofile->close();
    },
    outputs_definition => { step5_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step5_output file') },
    post_process_sub   => sub            { return 1 }
);

foreach my $step (@steps) {
    $multi_step_pipeline->add_step($step);
}
my (undef, undef, undef, $gid) = getpwuid $<;
my $default_group = getgrgid $gid;
my @groups = map { scalar getgrgid($_) } getgroups();
my $other_group;
foreach my $og (@groups) {
    next if $og eq $default_group;
    $other_group = $og;
    last;
}
my $second_pipelinesetup = VRPipe::PipelineSetup->create(name => 'ps2', datasource => $single_element_datasource, output_root => $second_pipeline_output_dir, pipeline => $multi_step_pipeline, $other_group ? (unix_group => $other_group) : ());

$single_element_datasource->elements;
my @second_output_files = (file(output_subdirs(VRPipe::DataElement->get(datasource => $single_element_datasource->id, result => { line => "single_element" })->id, 2), "1_step_1", 'output1.txt'));
foreach my $step_num (2 .. 5) {
    my $ofile = file(output_subdirs(VRPipe::DataElement->get(datasource => $single_element_datasource->id, result => { line => "single_element" })->id, 2), "${step_num}_step_$step_num", $step_num == 4 ? 'step4_basename.o' : "step${step_num}_output.txt");
    push(@second_output_files, $ofile);
}

# setup and use manager
ok my $manager = VRPipe::Manager->create(), 'can create a manager with no args';
is $manager->id, 1, 'manager always has an id of 1';
my @manager_setups = $manager->setups;
is @manager_setups, 2, 'setups() returns the correct number of PipelineSetups';

is handle_pipeline(@first_output_files, @second_output_files), 1, 'multi-step pipeline completed via Manager';

# check that the file group and directory permissions are as expected
SKIP: {
    skip "You do not belong to more than 1 group, so can't test changing group permissions", 3 unless $other_group;
    ok check_group($default_group, @first_output_files),  'When unix_group is not set, output files belong to default group';
    ok check_group($other_group,   @second_output_files), 'When unix_group is set, output files belong to the desired group';
    
    my $perms_correct = 1;
    foreach my $path (@second_output_files) {
        my $file = file($path);
        my $dir  = $file->dir;
        my (undef, undef, $mode) = stat($dir);
        $mode = sprintf("%04o", $mode & 07777);
        unless ($mode eq '0755') {
            $perms_correct = 0;
            last;
        }
    }
    ok $perms_correct, 'When unix_group is set, the directories get set to 0755';
    
    sub check_group {
        my (undef, undef, $desired_gid) = getgrnam(shift);
        foreach my $file (@_) {
            my (undef, undef, $mode, undef, undef, $gid) = stat($file);
            return 0 if $gid != $desired_gid;
        }
        return 1;
    }
}

# now lets create a pipeline using a pre-written step, where we'll test that a
# step can work with both a datasource input and the outputs of a previous step
my $prewritten_step_pipeline = VRPipe::Pipeline->create(name => 'md5_pipeline', description => 'simple test pipeline with two steps that create md5 files');
VRPipe::StepAdaptor->create(pipeline => $prewritten_step_pipeline, to_step => 1, adaptor_hash => { md5_file_input => { data_element => 0 } });
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));
VRPipe::StepAdaptor->create(pipeline => $prewritten_step_pipeline, to_step => 2, adaptor_hash => { md5_file_input => { md5_files => 1 } });
$prewritten_step_pipeline->add_step(VRPipe::Step->get(name => "md5_file_production"));

my $fofn_datasource = VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.fofn))->absolute);
$fofn_datasource->elements;
my $prewritten_step_pipeline_output_dir = dir($output_root, 'md5_pipeline');
my $md5_pipelinesetup = VRPipe::PipelineSetup->create(
    name        => 'ps3',
    datasource  => $fofn_datasource,
    output_root => $prewritten_step_pipeline_output_dir,
    pipeline    => $prewritten_step_pipeline,
    options     => { md5_files_in_source_dir => 0 }
);

my @md5_output_files = (file(output_subdirs(7, 3), '1_md5_file_production', 'file.bam.md5'), file(output_subdirs(8, 3), '1_md5_file_production', 'file.cat.md5'), file(output_subdirs(9, 3), '1_md5_file_production', 'file.txt.md5'), file(output_subdirs(7, 3), '2_md5_file_production', 'file.bam.md5.md5'), file(output_subdirs(8, 3), '2_md5_file_production', 'file.cat.md5.md5'), file(output_subdirs(9, 3), '2_md5_file_production', 'file.txt.md5.md5'));

is handle_pipeline(@md5_output_files), 1, 'all md5 files were created via Manager';

my @md5s;
foreach my $file (file(qw(t data file.bam))->absolute, file(qw(t data file.cat))->absolute, file(qw(t data file.txt))->absolute, file(output_subdirs(7, 3), '1_md5_file_production', 'file.bam.md5'), file(output_subdirs(8, 3), '1_md5_file_production', 'file.cat.md5'), file(output_subdirs(9, 3), '1_md5_file_production', 'file.txt.md5')) {
    push(@md5s, VRPipe::File->get(path => $file)->md5);
}
is_deeply [@md5s[0 .. 2]], [qw(21efc0b1cc21390f4dcc97795227cdf4 2f8545684149f81e26af90dec0c6869c eb8fa3ffb310ce9a18617210572168ec)], 'md5s were all set in db';

ok $md5s[3] && $md5s[4] && $md5s[5], 'md5s of md5 files were all set in db';

# test that a new pipelinesetup will work when using a previously completed
# datasource
$prewritten_step_pipeline_output_dir = dir($output_root, 'md5_pipeline2');
$md5_pipelinesetup = VRPipe::PipelineSetup->create(
    name        => 'ps4',
    datasource  => $fofn_datasource,
    output_root => $prewritten_step_pipeline_output_dir,
    pipeline    => $prewritten_step_pipeline,
    options     => { md5_files_in_source_dir => 0 }
);

@md5_output_files = (file(output_subdirs(7, 4), '1_md5_file_production', 'file.bam.md5'), file(output_subdirs(8, 4), '1_md5_file_production', 'file.cat.md5'), file(output_subdirs(9, 4), '1_md5_file_production', 'file.txt.md5'), file(output_subdirs(7, 4), '2_md5_file_production', 'file.bam.md5.md5'), file(output_subdirs(8, 4), '2_md5_file_production', 'file.cat.md5.md5'), file(output_subdirs(9, 4), '2_md5_file_production', 'file.txt.md5.md5'));

is handle_pipeline(@md5_output_files), 1, 'all md5 files were created via Manager, using a previously completed datasource';

# test that the step behaviour delete_inputs() won't delete input files needed
# by other setups that haven't finished yet
{
    my %initial_sofs;
    foreach my $f (@first_output_files) {
        my $file = VRPipe::File->get(path => $f);
        my ($sof) = VRPipe::StepOutputFile->search({ file => $file->id });
        $initial_sofs{ $file->id } = $sof->id;
    }
    
    my $sleep_step = VRPipe::Step->create(
        name              => 'sleep_step',
        inputs_definition => {},
        body_sub          => sub {
            my $self  = shift;
            my $sleep = int(rand(120)) + 60;
            $self->dispatch([qq{sleep $sleep}, $self->new_requirements(memory => 60, time => 1)]);
        },
        outputs_definition => {},
        post_process_sub   => sub { return 1 },
        description        => 'sleep for a variable amount of time'
    );
    
    my $step = VRPipe::Step->create(
        name              => 'file_input_to_output',
        inputs_definition => { file_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'file input', max_files => -1) },
        body_sub          => sub {
            my $self        = shift;
            my @input_files = @{ $self->inputs->{file_input} };
            my $basename    = join('-', map { $_->basename } @input_files);
            my $ofile       = $self->output_file(output_key => 'file_output', basename => "$basename.fito", type => 'txt')->path;
            $self->dispatch([qq{echo "output for $basename" > $ofile}, $self->new_requirements(memory => 60, time => 1)]);
        },
        outputs_definition => { file_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'fito file') },
        post_process_sub   => sub           { return 1 },
        description        => 'outputs a file corresponding to the input file'
    );
    
    my $fito_pipeline = VRPipe::Pipeline->create(name => 'fito_pipeline', description => 'file input to output pipeline');
    $fito_pipeline->add_step($sleep_step);
    $fito_pipeline->add_step($step);
    VRPipe::StepAdaptor->create(pipeline => $fito_pipeline, to_step => 2, adaptor_hash => { file_input => { data_element => 0 } });
    VRPipe::StepBehaviour->create(pipeline => $fito_pipeline, after_step => 2, behaviour => 'delete_inputs', behaviour_array => [['delete_inputs', 0]], default_regulation => 1);
    
    $step = VRPipe::Step->create(
        name              => 'file_input_to_output_2',
        inputs_definition => { file_input => VRPipe::StepIODefinition->create(type => 'txt', description => 'file input') },
        body_sub          => sub {
            my $self         = shift;
            my ($input_file) = @{ $self->inputs->{file_input} };
            my $basename     = $input_file->basename;
            my $ofile        = $self->output_file(output_key => 'file_output', basename => "$basename.fito2", type => 'txt')->path;
            $self->dispatch([qq{echo "output for $basename 2" > $ofile}, $self->new_requirements(memory => 60, time => 1)]);
        },
        outputs_definition => { file_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'fito2 file') },
        post_process_sub   => sub           { return 1 },
        description        => 'outputs a file corresponding to the input file'
    );
    
    my $fito2_pipeline = VRPipe::Pipeline->create(name => 'fito2_pipeline', description => 'file input to output pipeline 2');
    $fito2_pipeline->add_step($sleep_step);
    $fito2_pipeline->add_step($step);
    VRPipe::StepAdaptor->create(pipeline => $fito2_pipeline, to_step => 2, adaptor_hash => { file_input => { data_element => 0 } });
    VRPipe::StepBehaviour->create(pipeline => $fito2_pipeline, after_step => 2, behaviour => 'delete_inputs', behaviour_array => [['delete_inputs', 0]], default_regulation => 1);
    
    my $fito_output_dir = dir($output_root, 'fito_pipeline');
    
    my $vrds = VRPipe::DataSource->create(
        type   => 'vrpipe',
        method => 'all',
        source => 'ps1[1]'
    );
    my $fito_ps_1 = VRPipe::PipelineSetup->create(name => 'fito_ps1', datasource => $vrds, output_root => $fito_output_dir, pipeline => $fito_pipeline);
    my $fito_ps_2 = VRPipe::PipelineSetup->create(name => 'fito_ps2', datasource => $vrds, output_root => $fito_output_dir, pipeline => $fito2_pipeline);
    my $vrds2     = VRPipe::DataSource->create(
        type   => 'vrpipe',
        method => 'group_all',
        source => 'ps1[1]'
    );
    my $fito_ps_3 = VRPipe::PipelineSetup->create(name => 'fito_ps1', datasource => $vrds2, output_root => $fito_output_dir, pipeline => $fito_pipeline);
    
    my @fito_output_files;
    my $pager = $vrds->elements;
    while (my $elements = $pager->next) {
        my $i = 1;
        foreach my $el (@$elements) {
            push(@fito_output_files, file(output_subdirs($el->id, $fito_ps_1->id), '2_file_input_to_output',   "fed_result_$i.o.fito"));
            push(@fito_output_files, file(output_subdirs($el->id, $fito_ps_2->id), '2_file_input_to_output_2', "fed_result_$i.o.fito2"));
            $i++;
        }
    }
    $pager = $vrds2->elements;
    my ($el) = @{ $pager->next };
    push(@fito_output_files, file(output_subdirs($el->id, $fito_ps_3->id), '2_file_input_to_output', "fed_result_1.o-fed_result_2.o-fed_result_3.o-fed_result_4.o-fed_result_5.o.fito"));
    
    ok handle_pipeline(@fito_output_files), '3 setups that use the same input files and delete their inputs all ran successfully';
    my $ok = 0;
    foreach my $f (@first_output_files) {
        next if -e $f;
        my $file = VRPipe::File->get(path => $f);
        my ($sof) = VRPipe::StepOutputFile->search({ file => $file->id });
        if ($sof->id == $initial_sofs{ $file->id }) {
            $ok++;
        }
    }
    is $ok, 5, 'the input files got deleted, and were not recreated at any point';
}

# we can have the input step to a pipeline take 2 different types of file from
# the datasource
{
    my $ds = VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data datasource.fofnwm_mixed_types)),
        options => { metadata_keys => 'lane' }
    );
    
    # make a single-step pipeline that takes 2 different types as input
    my $single_step = VRPipe::Step->create(
        name              => 'two_type_step',
        inputs_definition => {
            bam_files   => VRPipe::StepIODefinition->create(type => 'bam', description => 'bam files',   max_files => -1, metadata => { lane => 'lane name' }),
            fastq_files => VRPipe::StepIODefinition->create(type => 'fq',  description => 'fastq files', max_files => -1, metadata => { lane => 'lane name' })
        },
        body_sub => sub {
            my $self = shift;
            my @bams = @{ $self->inputs->{bam_files} };
            my @fqs  = @{ $self->inputs->{fastq_files} };
            if (@bams == 1 && @fqs == 2) {
                foreach my $file (@bams, @fqs) {
                    $file->add_metadata({ ok => 1 });
                }
            }
        },
        outputs_definition => {},
        post_process_sub   => sub { return 1 },
        description        => 'a step that takes 2 different file types'
    );
    my $two_type_pipeline = VRPipe::Pipeline->create(name => 'two_type_step_pipeline', description => 'two_type_step pipeline');
    VRPipe::StepAdaptor->create(pipeline => $two_type_pipeline, to_step => 1, adaptor_hash => { bam_files => { data_element => 0 }, fastq_files => { data_element => 0 } });
    $two_type_pipeline->add_step($single_step);
    
    my $output_root = get_output_dir('manager_null_output');
    VRPipe::PipelineSetup->create(name => 'two_type_one_step_ps', datasource => $ds, output_root => $output_root, pipeline => $two_type_pipeline);
    handle_pipeline();
    
    my $oks = 0;
    foreach my $basename (qw(2822_6_1.fastq 2822_6_2.fastq 2822_6.pe.bam 2822_7_1.fastq 2822_7_2.fastq 2822_7.pe.bam)) {
        my $file = VRPipe::File->get(path => file('t', 'data', $basename)->absolute);
        $oks++ if $file->metadata->{ok};
    }
    is $oks, 6, 'we were able to run a pipeline where a step took files of 2 different types from the datasource';
    
    # it should also work with arbitrary file types where we have no
    # VRPipe::FileType::module written
    $ds = VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data datasource.fofnwm_mixed_unknown_types)),
        options => { metadata_keys => 'lane' }
    );
    $single_step = VRPipe::Step->create(
        name              => 'two_unknown_type_step',
        inputs_definition => {
            type1_files => VRPipe::StepIODefinition->create(type => 'typ1', description => 'typ1 files', max_files => -1, metadata => { lane => 'lane name' }),
            type2_files => VRPipe::StepIODefinition->create(type => 'typ2', description => 'typ2 files', max_files => -1, metadata => { lane => 'lane name' })
        },
        body_sub => sub {
            my $self   = shift;
            my @type1s = @{ $self->inputs->{type1_files} };
            my @type2s = @{ $self->inputs->{type2_files} };
            if (@type2s == 1 && @type1s == 2) {
                foreach my $file (@type1s, @type2s) {
                    $file->add_metadata({ ok => 1 });
                }
            }
        },
        outputs_definition => {},
        post_process_sub   => sub { return 1 },
        description        => 'a step that takes 2 different unknown file types'
    );
    $two_type_pipeline = VRPipe::Pipeline->create(name => 'two_unknown_type_step_pipeline', description => 'two_unknown_type_step pipeline');
    VRPipe::StepAdaptor->create(pipeline => $two_type_pipeline, to_step => 1, adaptor_hash => { type1_files => { data_element => 0 }, type2_files => { data_element => 0 } });
    $two_type_pipeline->add_step($single_step);
    
    VRPipe::PipelineSetup->create(name => 'two_unknown_type_one_step_ps', datasource => $ds, output_root => $output_root, pipeline => $two_type_pipeline);
    handle_pipeline();
    
    $oks = 0;
    foreach my $basename (qw(2822_6_1.typ1 2822_6_2.typ1 2822_6.pe.typ2 2822_7_1.typ1 2822_7_2.typ1 2822_7.pe.typ2)) {
        my $file = VRPipe::File->get(path => file('t', 'data', $basename)->absolute);
        $oks++ if $file->metadata->{ok};
    }
    is $oks, 6, 'we were able to run a pipeline where a step took files of 2 different unknown types from the datasource';
}

# when a job fails and is retried, test that we can get access to previous
# scheduler and job stdout/err
{
    my ($output_root, $pipeline) = create_single_step_pipeline('test_step_fail', 'file');
    my $fail_ps = VRPipe::PipelineSetup->create(name => 'fail_ps', datasource => $fofn_datasource, output_root => $output_root, pipeline => $pipeline);
    my $ps_id   = $fail_ps->id;
    my @ofiles  = (file(output_subdirs(7, $ps_id), '1_test_step_fail', 'file.bam'), file(output_subdirs(8, $ps_id), '1_test_step_fail', 'file.cat'), file(output_subdirs(9, $ps_id), '1_test_step_fail', 'file.txt'));
    
    ok handle_pipeline(@ofiles), 'pipeline with a step that fails twice before working ran successfully';
    
    my @subs = VRPipe::Submission->search({ 'stepstate.pipelinesetup' => $ps_id }, { order_by => { -asc => 'me.id' }, join => ['stepstate'] });
    my %job_std;
    sleep(5); # give some time for -handlers to call archive_output()
    foreach my $sub (@subs) {
        my $pars     = $sub->job_stdout;
        my $pr       = $pars->parsed_record;
        my $line_num = 0;
        while ($pars->next_record) {
            $line_num++;
            my $line = $pr->[0] || 'undef';
            $job_std{out}->{ $line_num . ' - ' . $line }++;
            last if $line_num == 3;
        }
        
        undef $pars;
        $pars     = $sub->job_stderr;
        $pr       = $pars->parsed_record;
        $line_num = 0;
        while ($pars->next_record) {
            $line_num++;
            my $line = $pr->[0] || 'undef';
            $job_std{err}->{ $line_num . ' - ' . $line }++;
            last if $line_num == 3;
        }
    }
    is_deeply \%job_std,
      {
        out => { "1 - undef" => 3, "2 - stdout message: failing on purpose since this is try 2" => 3, "3 - stdout message: failing on purpose since this is try 1" => 3 },
        err => { "1 - undef" => 3, "2 - stderr message: failing on purpose since this is try 2" => 3, "3 - stderr message: failing on purpose since this is try 1" => 3 }
      },
      'The job stdout and stderr of all 3 attempts on all 3 elements could be retrieved';
}

# test vrpipe-submissions --full_reset and --partial_reset work as expected
{
    # make a single-step pipeline that fails or succeeds based on metadata
    my $output_root = get_output_dir('fail_on_demand_pipeline');
    
    my $ds = VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data annotation.fofn))
    );
    
    my $input_file = VRPipe::File->create(path => file(qw(t data annotation.vcf.gz))->absolute, metadata => { fail => 1 });
    
    my $step = VRPipe::Step->create(
        name              => 'fail_on_demand',
        inputs_definition => { fod_input => VRPipe::StepIODefinition->create(type => 'vcf', description => 'fod input') },
        body_sub          => sub {
            my $self         = shift;
            my ($input_file) = @{ $self->inputs->{fod_input} };
            my $fail         = $input_file->meta_value('fail');
            my $req          = $self->new_requirements(memory => 1, time => 1);
            for (1 .. 3) {
                my $ofile = $self->output_file(output_key => 'fod_output', basename => "output$_.txt", type => 'txt');
                my $opath = $ofile->path;
                my $cmd   = qq{echo "output for $_" >> $opath};
                if ($_ == 2 && $fail) {
                    $self->dispatch([qq{$cmd; false}, $req, { output_files => [$ofile] }]);
                }
                else {
                    $self->dispatch([qq{$cmd; true}, $req]);
                }
            }
        },
        outputs_definition => { fod_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'fod_output file') },
        post_process_sub   => sub          { return 1 },
        description        => 'fail on demand step'
    );
    
    my $pipeline = VRPipe::Pipeline->create(name => 'fod_pipeline', description => 'pipeline that fails on demand');
    VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => 1, adaptor_hash => { fod_input => { data_element => 0 } });
    $pipeline->add_step($step);
    
    my $setup = VRPipe::PipelineSetup->create(name => 'fod setup', datasource => $ds, output_root => $output_root, pipeline => $pipeline);
    wait_till_setup_done($setup);
    
    my (@ofiles, @mtimes);
    my ($element) = VRPipe::DataElement->search({ datasource => $ds });
    foreach my $i (1 .. 3) {
        push(@ofiles, file(output_subdirs($element->id, $setup->id), '1_fail_on_demand', "output$i.txt"));
        push(@mtimes, $ofiles[-1]->stat->mtime);
    }
    
    my $perl     = VRPipe::Interface::CmdLine->vrpipe_perl_command('testing');
    my $setup_id = $setup->id;
    `$perl scripts/vrpipe-submissions --deployment testing --setup $setup_id -f --full_reset`;
    
    wait_till_setup_done($setup);
    
    my $all_good = 1;
    for (0 .. 2) {
        my $ofile = $ofiles[$_];
        
        my $expected_content = "output for " . ($_ + 1) . "\n";
        if ($expected_content ne $ofile->slurp) {
            $all_good = 0;
            last;
        }
        
        my $mtime = $ofile->stat->mtime;
        unless ($mtime > $mtimes[$_]) {
            $all_good = 0;
            last;
        }
        $mtimes[$_] = $mtime;
    }
    ok $all_good, 'vrpipe-submissions -f --full_reset worked as expected';
    
    `$perl scripts/vrpipe-submissions --deployment testing --setup $setup_id --partial_reset`;
    
    wait_till_setup_done($setup);
    
    $all_good = 1;
    for (0 .. 2) {
        my $ofile = $ofiles[$_];
        my $mtime = $ofile->stat->mtime;
        
        if ($_ == 1) {
            unless ($mtime > $mtimes[1]) {
                $all_good = 0;
                last;
            }
        }
        else {
            unless ($mtime == $mtimes[$_]) {
                $all_good = 0;
                last;
            }
        }
        
        my $expected_content = "output for " . ($_ + 1) . "\n";
        if ($expected_content ne $ofile->slurp) {
            $all_good = 0;
            last;
        }
    }
    ok $all_good, 'vrpipe-submissions --partial_reset worked as expected';
}

finish;
