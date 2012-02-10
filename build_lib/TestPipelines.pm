package TestPipelines;
use strict;
use warnings;
use Exporter 'import';
use Path::Class;
use lib "t";

our @EXPORT = qw(get_output_dir handle_pipeline output_subdirs create_single_step_pipeline get_bam_header get_bam_records finish);

our $manager = VRPipe::Manager->get();
our $scheduler;
our %setups;

BEGIN {
    $scheduler = VRPipe::Scheduler->get();
    $scheduler->start_scheduler;
}

sub get_output_dir {
    my $sub_dir = shift;
    my $output_dir = dir($scheduler->output_root, 'pipelines_test_output', $sub_dir);
    $scheduler->remove_tree($output_dir);
    $scheduler->make_path($output_dir);
    return $output_dir;
}

sub handle_pipeline {
    my $give_up = 1000;
    my $max_retries = VRPipeTest::max_retries();
    my $debug = VRPipeTest::debug();
    $manager->set_verbose_global(1) if $debug;
    while (1) {
        last if ($manager->trigger && all_pipelines_started());
        last if $give_up-- <= 0;
        $manager->handle_submissions(max_retries => $max_retries);
        
        # check for repeated failures
        my $submissions = $manager->unfinished_submissions();
        if ($submissions) {
            foreach my $sub (@$submissions) {
                next unless $sub->failed;
                if ($sub->retries >= $max_retries) {
                    warn "some submissions failed ", ($max_retries + 1), " times, giving up\n";
                    $manager->set_verbose_global(0) if $debug;
                    return 0;
                }
            }
        }
        
        sleep(1);
    }
    
    my $all_created = 1;
    foreach my $ofile (@_) {
        unless (-s $ofile) {
            warn "$ofile is missing\n";
            $all_created = 0;
        }
    }
    
    $manager->set_verbose_global(0) if $debug;
    return $all_created;
}

sub output_subdirs {
    my $de_id = shift;
    my $setup_id = shift || 1;
    my $setup = $setups{$setup_id};
    unless ($setup) {
        $setup = VRPipe::PipelineSetup->get(id => $setup_id);
        $setup->datasource->incomplete_element_states($setup); # create all dataelements et al.
        $setups{$setup_id} = $setup;
    }
    my $pipeline_root = $setup->output_root;
    
    my $des_id = VRPipe::DataElementState->get(dataelement => $de_id, pipelinesetup => $setup)->id;
    my $hashing_string = 'VRPipe::DataElementState::'.$des_id;
    my @subdirs = $manager->hashed_dirs($hashing_string);
    
    return ($pipeline_root, @subdirs, $de_id);
}

sub create_single_step_pipeline {
    my ($step_name, $input_key) = @_;
    
    my $step = VRPipe::Step->get(name => $step_name) || die "Could not create a step named '$step_name'\n";
    my $pipeline_name = $step_name.'_pipeline';
    my $pipeline = VRPipe::Pipeline->get(name => $pipeline_name, description => 'test pipeline for the '.$step_name.' step');
    VRPipe::StepMember->get(step => $step, pipeline => $pipeline, step_number => 1);
    VRPipe::StepAdaptor->get(pipeline => $pipeline, to_step => 1, adaptor_hash => { $input_key => { data_element => 0 } });
    
    my $output_dir = get_output_dir($pipeline_name);
    
    return ($output_dir, $pipeline, $step);
}

sub all_pipelines_started {
    my @setups = $manager->setups;
    my $schema = $manager->result_source->schema;
    foreach my $setup (@setups) {
        my $rs = $schema->resultset('DataElement')->search({ datasource => $setup->datasource->id, withdrawn => 0 });
        return 0 unless $rs->next;
    }
    return 1;
}

sub get_bam_header {
    my $bam = shift;
    my $bam_path = $bam->absolute;
    my $header = `samtools view -H $bam_path`;
    return split /\n/, $header;
}

sub get_bam_records {
    my $bam = shift;
    my $bam_path = $bam->absolute;
    my $records = `samtools view $bam_path`;
    return split /\n/, $records;
}

sub finish {
    $scheduler->stop_scheduler;
    exit;
}

1;