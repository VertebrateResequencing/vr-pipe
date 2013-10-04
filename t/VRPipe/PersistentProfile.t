#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use TryCatch;
use Time::HiRes qw(gettimeofday tv_interval);
use DateTime;

BEGIN {
    use Test::Most tests => 10;
    # only the author needs to run this test
    use VRPipeTest (required_env => 'VRPIPE_AUTHOR_TESTS');
    use TestPipelines;
    
    use_ok('VRPipe::Persistent::Schema');
    use_ok('VRPipe::Persistent::InMemory');
}

my %times;
my %elapsed;
my $l;
my $job_offset = 0;

# job creation tests
my $j1 = VRPipe::Job->create(cmd => 'foo');
my $j2 = VRPipe::Job->get(cmd => 'foo');
is $j2->id, $j1->id, 'create followed by get worked';

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    VRPipe::Job->create(cmd => qq[echo "job $i ] . 'n' x 1000 . qq[";]);
}
elapsed($l, __LINE__);
$job_offset = 1001;

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    VRPipe::Job->get(cmd => qq[echo "job $i ] . 'n' x 1000 . qq[";]);
}
elapsed($l, __LINE__);

# file creation & metadata tests
$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    my $file = VRPipe::File->create(path => "/my/file/path/$i");
    $file->add_metadata({
            foo  => 'bar',
            baz  => 'loman',
            cat  => 'dog',
            fish => 'turnip'
        }
    );
}
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    $file->add_metadata({
            foo => 'boo',
            rat => 'king' . $i
        }
    );
}
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
my $correct_meta = 0;
for my $i (1 .. 1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
    my $metadata = $file->metadata;
    $correct_meta++ if ($metadata->{rat} eq 'king' . $i && $metadata->{foo} eq 'boo' && $metadata->{baz} eq 'loman');
}
elapsed($l, __LINE__);
is $correct_meta, 1000, 'basic file metadata test worked';

$l = start_clock(__LINE__);
for my $i (1 .. 1000) {
    my $file = VRPipe::File->get(path => "/my/file/path/$i");
}
elapsed($l, __LINE__);

# later we need a datasource row in the db to make the dataelements, so we just
# create a list one and then withdraw its single element
my $ds = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.onelist))->absolute->stringify);
$ds->elements;
my @elements = VRPipe::DataElement->search({});
$elements[0]->withdrawn(1);
$elements[0]->update;
my $ds_id = $ds->id;

# make a pipeline and setup we'll use later
my $pipeline = VRPipe::Pipeline->create(name => 'two_step_pipeline', description => 'simple test pipeline with 2 steps');
my @steps;
$steps[0] = VRPipe::Step->create(
    name               => 'step_1',
    options_definition => { reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file to map against') },
    inputs_definition  => {},
    body_sub           => sub {
        my $self    = shift;
        my $options = $self->options;
        my $ref     = Path::Class::file($options->{reference_fasta});
        my $ofile   = $self->output_file(output_key => 'step1_output', output_dir => $ref->dir->stringify, basename => $ref->basename . '.step1out', type => 'txt');
        my $fh      = $ofile->openw();
        print $fh "step1output\n";
        $ofile->close();
        $self->dispatch(["$ref => " . $ofile->path, $self->new_requirements(memory => 3900, time => 1), { block_and_skip_if_ok => 1 }]);
    },
    outputs_definition => { step1_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step1_output file') },
    post_process_sub   => sub            { return 1 },
    description        => 'the first step'
);
$steps[1] = VRPipe::Step->create(
    name              => "step_2",
    inputs_definition => {
        from_datasource => VRPipe::StepIODefinition->create(type => 'txt', max_files   => 2, description => 'from datasource'),
        from_step1      => VRPipe::StepIODefinition->create(type => 'txt', description => 'from step_1')
    },
    body_sub => sub {
        my $self  = shift;
        my $ofile = $self->output_file(output_key => 'step2_output', basename => 'step2_output.txt', type => 'txt');
        my $fh    = $ofile->openw();
        print $fh "step2output\n";
        $ofile->close();
        my @fastqs = map { $_->path } @{ $self->inputs->{from_datasource} };
        my ($processed_ref) = map { $_->path } @{ $self->inputs->{from_step1} };
        $self->dispatch(["$processed_ref + (@fastqs) => " . $ofile->path, $self->new_requirements(memory => 3900, time => 1), { output_files => [$ofile] }]);
    },
    outputs_definition => { step2_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'step2_output file') },
    post_process_sub   => sub            { return 1 }
);
VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => 2, adaptor_hash => { from_datasource => { data_element => 0 }, from_step1 => { step1_output => 1 } });
foreach my $step (@steps) {
    $pipeline->add_step($step);
}
my $output_root = get_output_dir('profiling_output_dir');
my $ref         = file($output_root, 'ref.fa')->stringify;
my $setup       = VRPipe::PipelineSetup->create(name => 'ps1', datasource => $ds, output_root => $output_root, pipeline => $pipeline, active => 0, options => { reference_fasta => $ref }, controlling_farm => 'testing_farm');

# dataelement (and their associated files and metadata) creation tests
# (pretending to be like sequence_index datasource since it is metadata heavy
#  and has multiple paths per dataelement)
$l = start_clock(__LINE__);
my $lanes_hash;
for my $i (1 .. 50) {
    my $lane       = 'lane_' . $i;
    my $reads      = int(rand(100000)) + 500000;
    my $bases      = $reads * 108;
    my $library    = 'library_' . int($i / 5);
    my $sample_id  = int($i / 10);
    my $sample     = 'study_' . $sample_id;
    my $population = 'population_' . int($i / 25);
    for my $j (1 .. 2) {
        my $meta = {
            lane           => $lane,
            study          => 'study',
            study_name     => 'study_name',
            center_name    => 'SC',
            sample_id      => $sample_id,
            sample         => $sample,
            population     => $population,
            platform       => 'ILLUMINA',
            library        => $library,
            reads          => $reads,
            bases          => $bases,
            analysis_group => 'low_cov',
            paired         => 1,
        };
        my $fastq = file($output_root, "${lane}_$j.fq");
        my $vrfile = VRPipe::File->create(path => $fastq, type => 'fq');
        unless ($vrfile->s) {
            my $fh = $vrfile->openw();
            print $fh "input file\n";
            $vrfile->close();
            $vrfile->update_stats_from_disc;
        }
        $vrfile->add_metadata($meta, replace_data => 0);
        push(@{ $lanes_hash->{$lane}->{paths} }, $fastq);
    }
}
elapsed($l, __LINE__);

my ($most_recent_element_id) = VRPipe::DataElement->get_column_values('id', { datasource => $ds_id }, { order_by => { -desc => ['id'] }, rows => 1 });
my @element_args;
foreach my $lane (sort keys %$lanes_hash) {
    my $hash_ref = $lanes_hash->{$lane};
    my $result_hash = { paths => $hash_ref->{paths}, lane => $lane };
    push(@element_args, { datasource => $ds_id, result => $result_hash });
}
$l = start_clock(__LINE__);
VRPipe::DataElement->search_rs({ datasource => $ds_id })->update({ withdrawn => 1 });
elapsed($l, __LINE__);
$l = start_clock(__LINE__);
VRPipe::DataElement->bulk_create_or_update(map { $_->{withdrawn} = 0 unless defined $_->{withdrawn}; $_; } @element_args);
elapsed($l, __LINE__);
$l = start_clock(__LINE__);
my @des_args;
my $pager = VRPipe::DataElement->get_column_values_paged('id', { datasource => $ds_id, $most_recent_element_id ? (id => { '>' => $most_recent_element_id }) : () });

while (my $eids = $pager->next) {
    foreach my $eid (@$eids) {
        push(@des_args, { pipelinesetup => $setup->id, dataelement => $eid });
    }
}
VRPipe::DataElementState->bulk_create_or_update(@des_args) if @des_args;
elapsed($l, __LINE__);

my @element_ids = VRPipe::DataElement->get_column_values('id', { withdrawn => 0 });
is_deeply \@element_ids, [2 .. 51], 'created correct dataelements';

# now test the speed of triggering all these dataelements through a 2 step
# pipeline, the first of which is a block_and_skip
$l = start_clock(__LINE__);
$setup->trigger(first_step_only => 1, prepare_elements => 0);
elapsed($l, __LINE__);

my @jobs = VRPipe::Job->search({ id => { '>' => $job_offset } });
my @subs = VRPipe::Submission->search({});
is_deeply [scalar(@jobs), scalar(@subs)], [1, 1], 'triggering first step only created 1 job and 1 sub';

# pretend we ran the submission successfully
sub pretend_completion {
    my ($jobs, $subs) = @_;
    foreach my $job (@$jobs) {
        $job->exit_code(0);
        foreach my $col (qw(start_time end_time)) {
            $job->$col(DateTime->now);
        }
        $job->update;
    }
    foreach my $sub (@$subs) {
        $sub->_done(1);
        $sub->update;
    }
}
pretend_completion(\@jobs, \@subs);

# trigger the second step
$l = start_clock(__LINE__);
unless ($setup->currently_complete) {
    my $ie_pager = $setup->datasource->incomplete_element_states($setup, prepare => 0);
    if ($ie_pager) {
        while (my $dess = $ie_pager->next) {
            foreach my $des (@$dess) {
                my $de = $des->dataelement;
                my $j  = start_clock(__LINE__);
                $setup->trigger(dataelement => $de);
                elapsed($j, __LINE__);
            }
        }
    }
}
elapsed($l, __LINE__);

@jobs = VRPipe::Job->search({ id => { '>' => $job_offset } });
@subs = VRPipe::Submission->search({});
is_deeply [scalar(@jobs), scalar(@subs)], [51, 51], 'forcing trigger on all dataelements created all jobs and subs';
$l = start_clock(__LINE__);
pretend_completion(\@jobs, \@subs);
elapsed($l, __LINE__);

# now withdraw some dataelements and add new ones
$l = start_clock(__LINE__);
undef $lanes_hash;
for my $i (6 .. 75) {
    my $lane       = 'lane_' . $i;
    my $reads      = int(rand(100000)) + 500000;
    my $bases      = $reads * 108;
    my $library    = 'library_' . int($i / 5);
    my $sample_id  = int($i / 10);
    my $sample     = 'study_' . $sample_id;
    my $population = 'population_' . int($i / 25);
    for my $j (1 .. 2) {
        my $meta = {
            lane           => $lane,
            study          => 'study',
            study_name     => 'study_name',
            center_name    => 'SC',
            sample_id      => $sample_id,
            sample         => $sample,
            population     => $population,
            platform       => 'ILLUMINA',
            library        => $library,
            reads          => $reads,
            bases          => $bases,
            analysis_group => 'low_cov',
            paired         => 1,
        };
        my $fastq = file($output_root, "${lane}_$j.fq");
        my $vrfile = VRPipe::File->create(path => $fastq, type => 'fq');
        unless ($vrfile->s) {
            my $fh = $vrfile->openw();
            print $fh "input file\n";
            $vrfile->close();
            $vrfile->update_stats_from_disc;
        }
        $vrfile->add_metadata($meta, replace_data => 0);
        push(@{ $lanes_hash->{$lane}->{paths} }, $fastq);
    }
}
elapsed($l, __LINE__);

($most_recent_element_id) = VRPipe::DataElement->get_column_values('id', { datasource => $ds_id }, { order_by => { -desc => ['id'] }, rows => 1 });
@element_args = ();
foreach my $lane (sort keys %$lanes_hash) {
    my $hash_ref = $lanes_hash->{$lane};
    my $result_hash = { paths => $hash_ref->{paths}, lane => $lane };
    push(@element_args, { datasource => $ds_id, result => $result_hash });
}
$l = start_clock(__LINE__);
VRPipe::DataElement->search_rs({ datasource => $ds_id })->update({ withdrawn => 1 });
elapsed($l, __LINE__);
$l = start_clock(__LINE__);
VRPipe::DataElement->bulk_create_or_update(map { $_->{withdrawn} = 0 unless defined $_->{withdrawn}; $_; } @element_args);
elapsed($l, __LINE__);
$l        = start_clock(__LINE__);
@des_args = ();
$pager    = VRPipe::DataElement->get_column_values_paged('id', { datasource => $ds_id, $most_recent_element_id ? (id => { '>' => $most_recent_element_id }) : () });

while (my $eids = $pager->next) {
    foreach my $eid (@$eids) {
        push(@des_args, { pipelinesetup => $setup->id, dataelement => $eid });
    }
}
VRPipe::DataElementState->bulk_create_or_update(@des_args) if @des_args;
elapsed($l, __LINE__);

my %element_ids = map { $_ => 1 } VRPipe::DataElement->get_column_values('id', { withdrawn => 0 });
is_deeply [scalar(keys %element_ids), exists $element_ids{75}], [70, 1], 'updated dataelements correctly';

# now retest the trigger, especially so we can test the speed of skipping things
# that have already been done
$l = start_clock(__LINE__);
$setup->trigger(first_step_only => 1, prepare_elements => 0);
elapsed($l, __LINE__);

$l = start_clock(__LINE__);
unless ($setup->currently_complete) {
    my $ie_pager = $setup->datasource->incomplete_element_states($setup, prepare => 0);
    if ($ie_pager) {
        while (my $dess = $ie_pager->next) {
            foreach my $des (@$dess) {
                my $de = $des->dataelement;
                my $j  = start_clock(__LINE__);
                $setup->trigger(dataelement => $de);
                elapsed($j, __LINE__);
            }
        }
    }
}
elapsed($l, __LINE__);

@jobs = VRPipe::Job->search({ id => { '>' => $job_offset } });
@subs = VRPipe::Submission->search({ 'dataelement.withdrawn' => 0 }, { join => { 'stepstate' => 'dataelement' } });
my %jobs = map { $_->job->id => 1 } @subs;
is_deeply [scalar(keys %jobs), scalar(@subs)], [70, 70], 'forcing trigger on all dataelements created correct jobs and subs';

# do what vrpipe-server does to work out what submissions need to be submitted
my $scheduler = VRPipe::Scheduler->create;
my $sched_id  = $scheduler->id;
$l = start_clock(__LINE__);
my $sub_pager = VRPipe::Submission->search_paged({ '_done' => 0, -or => [-and => ['_failed' => 1, retries => { '<' => 3 }], '_failed' => 0], scheduler => $sched_id, 'pipelinesetup.controlling_farm' => 'testing_farm', 'pipelinesetup.active' => 0, 'dataelement.withdrawn' => 0 }, { join => { stepstate => ['pipelinesetup', 'dataelement'] }, prefetch => ['requirements', 'job', { stepstate => { stepmember => 'step' } }] });

my $queued = 0;
my $im     = VRPipe::Persistent::InMemory->new();
while (my $subs = $sub_pager->next(no_resetting => 1)) {
    foreach my $sub (@$subs) {
        my $sub_id = $sub->id;
        # black_and_skip handling
        my $job = $sub->job;
        if ($job->block_and_skip_if_ok) {
            my ($first_sub) = VRPipe::Submission->search({ 'job' => $job->id, '_done' => 0 }, { rows => 1, order_by => { -asc => 'id' } });
            next unless ($first_sub && $first_sub->id == $sub_id);
        }
        
        #*** we're not testing global step limit handling, which might be slow
        
        # queue this one in redis
        my $req_id = $sub->requirements->id;
        $im->enqueue($req_id, $sub_id);
        $queued++;
    }
}
elapsed($l, __LINE__);
is $queued, 25, 'correctly queued 25 subs to be submitted';

report();

done_testing;
exit;

sub start_clock {
    my $l1 = shift;
    $times{$l1} = [gettimeofday];
    return $l1;
}

sub elapsed {
    my ($l1, $l2, $silent) = @_;
    my $e         = tv_interval($times{$l1});
    my $e_rounded = sprintf("%0.2f", $e);
    my $id        = "$l1..$l2";
    note("Going from line $id took $e_rounded seconds\n") unless $silent;
    push(@{ $elapsed{$id} }, $e);
}

sub report {
    my %cums;
    while (my ($id, $times) = each %elapsed) {
        my $cum = 0;
        foreach my $t (@$times) {
            $cum += $t;
        }
        $cums{$id} = $cum;
    }
    
    note("\nMost time consuming sections:\n");
    foreach my $id (sort { $cums{$b} <=> $cums{$a} || $a cmp $b } keys %cums) {
        my $cum = $cums{$id};
        
        my $times = $elapsed{$id};
        my $total = 0;
        my $count = 0;
        foreach my $t (@$times) {
            $count++;
            $total += $t;
        }
        
        my $avg = sprintf("%0.4f", $total / $count);
        $cum = sprintf("%0.2f", $cum);
        
        my $note = "\t$id: $cum seconds";
        if ($count > 1) {
            $note .= " ($avg avg over $count loops)\n";
        }
        else {
            $note .= "\n";
        }
        note($note);
    }
}

# MySQL, VRPipe v0.155 (introduction of keyvallists on DataElement and File):
# Most time consuming sections:
#   243..257: 167.70 seconds
#   250..252: 167.65 seconds (3.3529 avg over 50 loops)
#   341..355: 100.28 seconds
#   348..350: 100.20 seconds (1.0547 avg over 95 loops)
#   337..339: 83.83 seconds
#   56..65: 32.62 seconds [2x slower...]
#   43..54: 19.22 seconds [similar speed]
#   29..33: 13.00 seconds
#   36..40: 8.85 seconds
#   144..182: 7.13 seconds
#   67..74: 7.12 seconds [similar speed]
#   217..219: 6.35 seconds
#   267..305: 4.72 seconds
#   77..81: 3.84 seconds
#   194..196: 1.82 seconds
#   317..319: 1.59 seconds
#   262..264: 0.45 seconds
#   365..387: 0.08 seconds
#   320..330: 0.05 seconds
#   197..207: 0.04 seconds
#   314..316: 0.00 seconds
#   191..193: 0.00 seconds

# MySQL 5.5, VRPipe v0.164 (fixed adaptor_hash definition to avoid transaction
# failures):
# Most time consuming sections:
#   56..65: 28.28 seconds
#   29..33: 16.48 seconds
#   43..54: 15.58 seconds
#   36..40: 13.45 seconds
#   242..256: 12.56 seconds
#   249..251: 12.50 seconds (0.2501 avg over 50 loops)
#   336..338: 10.80 seconds
#   216..218: 9.84 seconds
#   340..354: 7.74 seconds
#   347..349: 7.66 seconds (0.0806 avg over 95 loops)
#   67..74: 6.67 seconds
#   143..181: 6.66 seconds
#   266..304: 4.94 seconds
#   77..81: 3.64 seconds
#   316..318: 1.87 seconds
#   193..195: 1.64 seconds
#   261..263: 0.28 seconds
#   364..386: 0.08 seconds
#   319..329: 0.04 seconds
#   196..206: 0.03 seconds
#   313..315: 0.00 seconds
#   190..192: 0.00 seconds

# MySQL 5.5, VRPipe v0.168 (after greater redis usage and connection fix)
# Most time consuming sections:
#   56..65: 31.40 seconds
#   239..253: 20.21 seconds
#   246..248: 20.16 seconds (0.4031 avg over 50 loops)
#   43..54: 19.33 seconds
#   333..335: 13.73 seconds
#   29..33: 12.70 seconds
#   36..40: 11.68 seconds
#   337..351: 11.43 seconds
#   344..346: 11.35 seconds (0.1195 avg over 95 loops)
#   213..215: 10.61 seconds
#   143..181: 7.05 seconds
#   67..74: 6.64 seconds
#   263..301: 5.07 seconds
#   77..81: 3.53 seconds
#   313..315: 1.85 seconds
#   193..195: 1.80 seconds
#   258..260: 0.24 seconds
#   361..384: 0.09 seconds
#   316..326: 0.04 seconds
#   196..206: 0.03 seconds
#   310..312: 0.00 seconds
#   190..192: 0.00 seconds
