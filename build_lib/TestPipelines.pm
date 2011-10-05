package TestPipelines;
use strict;
use warnings;
use TestPersistentReal;
use Exporter 'import';
use Path::Class;
use lib "t";

our @EXPORT = qw(get_output_dir handle_pipeline output_subdirs create_single_step_pipeline finish);

our $manager = VRPipe::Manager->get();
our $scheduler;

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
    while (! $manager->trigger) {
        last if $give_up-- <= 0;
        $manager->handle_submissions;
        
        # check for repeated failures
        my $submissions = $manager->unfinished_submissions();
        if ($submissions) {
            foreach my $sub (@$submissions) {
                next unless $sub->failed;
                if ($sub->retries >= 3) {
                    warn "some submissions failed 3 times, giving up\n";
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
    
    return $all_created;
}

sub output_subdirs {
    my $element_id = shift;
    return ($manager->hashed_dirs('VRPipe::DataElement::'.$element_id), $element_id);
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

sub finish {
    $scheduler->stop_scheduler;
    exit;
}

1;