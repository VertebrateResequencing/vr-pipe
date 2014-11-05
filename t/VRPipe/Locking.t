#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class qw(file dir);
use Parallel::ForkManager;
use Sys::Hostname;

BEGIN {
    use Test::Most tests => 9;
    use VRPipeTest;
    use TestPipelines;
}

# make a stepstate for some basic testing
my $test_pipeline = VRPipe::Pipeline->create(name => 'test_pipeline');
my $ds = VRPipe::DataSource->create(type => 'list', method => 'all', source => file(qw(t data datasource.onelist))->absolute);
my $ps = VRPipe::PipelineSetup->create(name => 'ps', datasource => $ds, output_root => dir(qw(tmp)), pipeline => $test_pipeline, options => {}, active => 0);
$ds->elements;
my $ss = VRPipe::StepState->create(stepmember => 1, dataelement => 1, pipelinesetup => 1);
my $des = VRPipe::DataElementState->get(id => 1);

# first just check that the basics of how we expect transactions and row locking
# to work
my $fm                 = Parallel::ForkManager->new(3);
my $did_complete       = 0;
my $already_complete   = 0;
my $confirmed_complete = 0;
$fm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my ($dc, $ac, $cc) = @$data_structure_reference;
        $did_complete++       if $dc;
        $already_complete++   if $ac;
        $confirmed_complete++ if $cc;
    }
);
for (1 .. 1) {
    $did_complete       = 0;
    $already_complete   = 0;
    $confirmed_complete = 0;
    $ss->reselect_values_from_db;
    $ss->complete(0);
    $ss->update;
    $des->reselect_values_from_db;
    $des->completed_steps(0);
    $des->update;
    
    foreach my $i (1 .. 3) {
        $fm->start and next;
        
        #sleep($i);
        
        $ss->block_until_locked;
        $ss->maintain_lock;
        
        my ($dc, $ac, $cc);
        #warn "$$ loop start, ss->complete is ", $ss->complete, "\n";
        my $transaction = sub {
            #warn "$$ transaction start\n";
            #my ($state) = VRPipe::StepState->search({ id => $ss->id }, { for => 'update' });
            #$state->reselect_values_from_db;
            
            if ($ss->complete) {
                #warn "$$ already complete\n";
                $ac = 1;
            }
            else {
                sleep(5);
                #warn "$$ completing\n";
                if (nested_transaction()) {
                    #warn "$$ nested_transaction returned true\n";
                    $des->reselect_values_from_db;
                    if ($des->completed_steps == 1) {
                        #warn "$$ des completed_steps was 1\n";
                        $ss->complete(1);
                        $ss->update;
                        $dc = 1;
                        #warn "$$ set ss->complete(1)\n";
                        
                        if (nested_transaction_two()) {
                            #warn "$$ nested_transaction_two also returned true, will set cc = 1\n";
                            $cc = 1;
                        }
                    }
                    else {
                        #warn "$$ des->completed_steps was ", $des->completed_steps, "\n";
                    }
                }
            }
            #warn "$$ transaction end\n";
        };
        $ss->do_transaction($transaction, 'failed');
        $ss->unlock;
        #warn "$$ loop end, ss->complete is ", $ss->complete, "\n";
        $fm->finish(0, [$dc, $ac, $cc]);
    }
    $fm->wait_all_children;
    is_deeply [$did_complete, $already_complete, $confirmed_complete], [1, 2, 1], 'basics of transactions and row locking work';
}

sub nested_transaction {
    my $at_zero = 0;
    
    my $des = VRPipe::DataElementState->get(id => 1);
    #warn "   $$ will bul for des\n";
    $des->block_until_locked;
    #warn "   $$ got lock for des\n";
    
    my $transaction = sub {
        if ($des->completed_steps == 0) {
            $at_zero = 1;
            $des->completed_steps(1);
            $des->update;
        }
    };
    $ss->do_transaction($transaction, 'failed nested');
    $des->unlock;
    #warn "   $$ unlocked des, at_zero is $at_zero\n";
    return $at_zero;
}

sub nested_transaction_two {
    my $found       = 0;
    my $transaction = sub {
        ($found) = VRPipe::StepState->get_column_values('id', { complete => 1 });
    };
    $ss->do_transaction($transaction, 'failed nested 2');
    return $found;
}

# now do a more real-world test that used to reveal an issue that broke the
# above locking assumptions and allowed multiple processes to work on the same
# stepstate at the same time, causing random havoc
my ($output_dir, $pipeline, $step) = create_single_step_pipeline('fake_fastq_metadata', 'fastq_files');
my $si_datasource = VRPipe::DataSource->create(
    type    => 'sequence_index',
    method  => 'lane_fastqs',
    source  => file(qw(t data datasource.sequence_index))->absolute,
    options => { local_root_dir => dir(".")->absolute->stringify }
);
my $setup = VRPipe::PipelineSetup->create(
    name        => 'fm_setup',
    datasource  => $si_datasource,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);

# make 2 more setups in different dirs; whilst we don't have any tests
# specific to these setups, creating them is sufficient to reveal problems
my $output_dir2 = get_output_dir('fm_setup2');
VRPipe::PipelineSetup->create(
    name        => 'fm_setup2',
    datasource  => $si_datasource,
    output_root => $output_dir2,
    pipeline    => $pipeline,
    options     => {}
);
my $output_dir3 = get_output_dir('fm_setup3');
VRPipe::PipelineSetup->create(
    name        => 'fm_setup3',
    datasource  => $si_datasource,
    output_root => $output_dir3,
    pipeline    => $pipeline,
    options     => {}
);

my @ofiles;
foreach my $basename (qw(2822_6.fastq 2822_6_1.fastq 2822_6_2.fastq 2822_7_1.fastq 2822_7_2.fastq 2823_4_1.fastq 2823_4_2.fastq 8324_8_1.fastq 8324_8_2.fastq)) {
    push(@ofiles, file('t', 'data', $basename)->absolute);
}

ok handle_pipeline(@ofiles), 'fastq_metadata pipeline ran ok and all input files still exist';

my $fastq = VRPipe::File->create(path => file(qw(t data 2822_6.fastq))->absolute);
is $fastq->metadata->{avg_read_length}, 10, 'it actually added some metadata';

# the main symptom of problems is that start_over happens
my @logs = VRPipe::PipelineSetupLog->search({ message => { like => "%StepState->start_over was called%" } });
is scalar(@logs), 0, 'no start_over calls were made';

# start_over, and other things going wrong, happen due to Jobs running more than
# once at a time, or triggers on the same stepstate happening more than once
# at a time. we don't mind the same pid repeating these things, since that
# suggests just a retry of a loop (serial repetition) instead of something
# actually bad happening (parallel repetition)
check_jobs_ran_once();

my @ss_ids = VRPipe::StepState->get_column_values('id', { pipelinesetup => { '!=' => 1 } });
my $triggered_once_count = 0;
foreach my $ss_id (@ss_ids) {
    my @logs = VRPipe::PipelineSetupLog->search({ ss_id => $ss_id, message => { like => "%so will trigger the next Step%" } });
    my %by_pid = map { $_->pid => 1 } @logs;
    my $count = keys %by_pid;
    $triggered_once_count++ if $count <= 1; # this should be == 1, but db issues mean we don't always create the Log rows we requested
    if ($count > 1) {                       # this should be != 1
        warn "ss_id $ss_id is bad\n";
    }
}
is $triggered_once_count, scalar(@ss_ids), 'all stepstates only triggered the next step once';

# following test for author only on LSF, since we need 1000+ submissions running
# at once
SKIP: {
    my $host = hostname();
    skip "author-only 1000 submission tests", 3 unless $host eq 'vr-2-2-02' && $ENV{VRPIPE_OPTIONAL_TESTS};
    
    # long after the above problems were resolved, it was found that a bam_index
    # pipeline fell over when given a datasource of ~3000 dataelements where
    # there were groups of ~50 dataelements that all had the same list of ~1000
    # bams. The symptoms were that tons of select for updates on the same
    # stepstate id would clog up all the database connections, and the same job
    # would run multiple times in a row, exiting 0 but being considered dead by
    # the next handler to claim_and_run it
    
    # confirmed that this test at least duplicated the mass (>600) of processes
    # all doing:
    # SELECT me.pipelinesetup, me.dataelement, me.id, me.cmd_summary, me.same_submissions_as, me.stepmember, me.complete FROM stepstate me WHERE ( id = '1' ) ORDER BY me.id ASC FOR UPDATE
    # this is while stepstates are slowly being created and we only had the
    # first 10 or so, all same_submissions_as 1 (and 1000 submissions)
    # the number of processes stacked up was similar to number of jobs running
    # in LSF. And it did soon stack over 1000, run out of connections and then
    # the test script died
    
    # fixing that issue then raised the problem of it being really slow (~10s+)
    # to create each stepstate, with the time spent mostly in running the
    # body_sub - for further improvement we need to farm out step handlers
    # instead of having the setup handler parse steps in a single-process loop
    
    # make a simple 1-step pipeline
    my $fake_bam_index_step = VRPipe::Step->create(
        name              => 'fake_bam_index',
        inputs_definition => { in_files => VRPipe::StepIODefinition->create(type => 'txt', description => 'input file', max_files => -1) },
        body_sub          => sub {
            my $self = shift;
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $in (@{ $self->inputs->{in_files} }) {
                my $in_path  = $in->path;
                my $out      = $self->output_file(output_key => 'out_files', output_dir => $in->dir, basename => $in->basename . '.out', type => 'txt');
                my $out_path = $out->path;
                $self->dispatch([qq{cat $in_path >> $out_path}, $req, { output_files => [$out], block_and_skip_if_ok => 1 }]);
            }
        },
        outputs_definition => { out_files => VRPipe::StepIODefinition->create(type => 'txt', description => 'output file', max_files => -1) },
        post_process_sub   => sub         { return 1 },
        description        => 'cat input to output'
    );
    
    my $pipeline = VRPipe::Pipeline->create(name => 'fake_bam_index_pipeline', description => 'simple test pipeline that behaves like bam_index');
    $pipeline->add_step($fake_bam_index_step);
    VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => 1, adaptor_hash => { in_files => { data_element => 0 } });
    
    # make the groups of input files we need to test with, and create
    # dataelements as we go since it is too slow to just call $datasource->
    # elements() afterwards in the normal way
    my $output_root = get_output_dir('locking_test');
    my $inputs_dir = dir($output_root, 'inputs');
    mkdir($inputs_dir);
    my $source_file = file($output_root, 'inputs.fofn');
    my $datasource = VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $source_file->stringify,
        options => { metadata_keys => "group" }
    );
    my $ds_id = $datasource->id;
    
    note "will make input files and dataelements...";
    my $source_fh = $source_file->openw;
    print $source_fh "path\tgroup\n";
    my @ofiles;
    my @e_args;
    for my $sub_dir (1 .. 3) {
        my $dir = dir($inputs_dir, $sub_dir);
        mkdir($dir);
        my $paths;
        for my $i (1 .. 1000) {
            my $ifile = file($dir, "file.$i");
            push(@$paths, $ifile->stringify);
            
            my $in_fh = $ifile->openw;
            print $in_fh "input file $i\n";
            close($in_fh);
            VRPipe::File->create(path => $ifile, type => 'txt');
            
            my $ofile = file($ifile . '.out');
            push(@ofiles, $ofile);
            
            for my $j (1 .. 10) {
                print $source_fh $ifile, "\t", "$sub_dir.$j", "\n";
            }
        }
        
        for my $j (1 .. 10) {
            push(@e_args, { datasource => $ds_id, result => { paths => $paths, group => "$sub_dir.$j" } });
        }
    }
    close($source_fh);
    VRPipe::DataElement->bulk_create_or_update(@e_args);
    note "... created dataelements, will run pipeline";
    
    my $sf = VRPipe::File->create(path => $source_file);
    $sf->update_md5;
    $datasource->_changed_marker($sf->md5);
    $datasource->update;
    $datasource->elements;
    
    my $pipelinesetup = VRPipe::PipelineSetup->create(
        name        => 'ps1',
        datasource  => $datasource,
        output_root => $output_root,
        pipeline    => $pipeline
    );
    
    # create and run the setup and confirm all jobs only ran once
    ok handle_pipeline(@ofiles), 'fastq_metadata pipeline ran ok and all input files still exist';
    check_jobs_ran_once();
    my $out_lines = 0;
    foreach my $ofile (@ofiles) {
        my $fh    = $ofile->openr;
        my @lines = <$fh>;
        if (@lines != 1) {
            warn " $ofile had ", scalar(@lines), " lines\n";
        }
        $out_lines += @lines;
    }
    is $out_lines, 3000, 'each output file was only written to once';
}

done_testing;
exit;

sub check_jobs_ran_once {
    my @job_ids = VRPipe::Job->get_column_values('id', {});
    my $ran_once_count = 0;
    foreach my $jid (@job_ids) {
        my @logs = VRPipe::PipelineSetupLog->search({ job_id => $jid, message => { like => "%cmd-running child exited with code 0%" } });
        my %by_pid = map { $_->pid => 1 } @logs;
        my $count = keys %by_pid;
        $ran_once_count++ if $count == 1;
        if ($count != 1) {
            warn "job_id $jid is bad\n";
        }
    }
    is $ran_once_count, scalar(@job_ids), 'all jobs only ran once';
}
