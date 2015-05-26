
=head1 NAME

TestPipelines - test methods for use when testing pipelines

=head1 SYNOPSIS
    
    BEGIN {
        use Test::Most tests => 3;
        use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
                        required_exe => [qw(samtools bwa)]);
        use TestPipelines;
    }
    
    my $output_dir = get_output_dir('my_test_or_pipeline_name');
    
    # ... setup a pipeline ...
    
    ok handle_pipeline(@expected_files), 'pipeline ran ok and generated files';

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

package TestPipelines;
use strict;
use warnings;
use Exporter 'import';
use Path::Class;

our @EXPORT = qw(get_output_dir handle_pipeline wait_till_setup_done output_subdirs create_single_step_pipeline get_bam_header get_bam_records finish);

our $manager = VRPipe::Manager->create();
our $scheduler;
our %setups;

BEGIN {
    $scheduler = VRPipe::Scheduler->create();
}

sub get_output_dir {
    my $sub_dir = shift;
    my $output_dir = dir($scheduler->output_root, 'pipelines_test_output', $sub_dir);
    $scheduler->remove_tree($output_dir);
    $scheduler->make_path($output_dir);
    return $output_dir;
}

sub handle_pipeline {
    my $give_up     = 3600;
    my $max_retries = VRPipeTest::max_retries();
    my $debug       = VRPipeTest::debug();
    $manager->set_verbose_global(1) if $debug;
    my $gave_up    = 0;
    my $retriggers = 0;
    while (1) {
        # check for repeated failures
        my $num_failed = VRPipe::Submission->search({ _failed => 1, retries => { '>=' => $max_retries } });
        if ($num_failed) {
            warn "$num_failed submissions failed ", ($max_retries + 1), " times, giving up\n";
            $manager->set_verbose_global(0) if $debug;
            return 0;
        }
        
        if (all_pipelines_finished()) {
            # make sure linked pipelines have a chance to get all their data
            # elements once their parent pipelines have completed
            $retriggers++;
            last if $retriggers >= 3;
        }
        
        if ($give_up-- <= 0) {
            $gave_up = 1;
            warn "not all pipelinesetups finished yet, but giving up after 1500 cycles\n";
            last;
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
    return $gave_up ? 0 : $all_created;
}

sub output_subdirs {
    my $de_id    = shift;
    my $setup_id = shift || 1;
    my $setup    = $setups{$setup_id};
    unless ($setup) {
        $setup = VRPipe::PipelineSetup->get(id => $setup_id);
        $setup->datasource->incomplete_element_states($setup); # create all dataelements et al.
        $setups{$setup_id} = $setup;
    }
    my $pipeline_root = $setup->output_root;
    
    my $des_id         = VRPipe::DataElementState->get(dataelement => $de_id, pipelinesetup => $setup)->id;
    my $hashing_string = 'VRPipe::DataElementState::' . $des_id;
    my @subdirs        = $manager->hashed_dirs($hashing_string);
    
    return ($pipeline_root, @subdirs, $de_id);
}

sub create_single_step_pipeline {
    my ($step_name, @input_keys) = @_;
    
    my $step          = VRPipe::Step->get(name => $step_name) || die "Could not create a step named '$step_name'\n";
    my $pipeline_name = $step_name . '_pipeline';
    my $pipeline      = VRPipe::Pipeline->create(name => $pipeline_name, description => 'test pipeline for the ' . $step_name . ' step');
    VRPipe::StepMember->create(step => $step, pipeline => $pipeline, step_number => 1);
    my %adapter_hash = map { $_ => { data_element => 0 } } @input_keys;
    VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => 1, adaptor_hash => \%adapter_hash);
    
    my $output_dir = get_output_dir($pipeline_name);
    
    return ($output_dir, $pipeline, $step);
}

sub all_pipelines_started {
    my @setups = $manager->setups;
    foreach my $setup (grep { $_->active } @setups) {
        # a test might change a datasource then immediately test for the
        # results, before the server calls elements in its watcher, so we call
        # it manually ourselves every time
        $setup->datasource->elements;
        
        my $found = VRPipe::DataElement->search({ datasource => $setup->datasource->id, withdrawn => 0 });
        return 0 unless $found;
    }
    return 1;
}

sub all_pipelines_finished {
    return 0 unless all_pipelines_started();
    
    my @setups = $manager->setups;
    foreach my $setup (grep { $_->active } @setups) {
        my $setup_id  = $setup->id;
        my $pipeline  = $setup->pipeline;
        my $num_steps = $pipeline->step_members;
        
        my $not_done = VRPipe::DataElementState->search({ pipelinesetup => $setup_id, 'dataelement.withdrawn' => 0, completed_steps => { '!=' => $num_steps } }, { join => 'dataelement' });
        
        return 0 if $not_done;
    }
    return 1;
}

sub wait_till_setup_done {
    # return true when the setup is done doing stuff - either completed
    # successfully or stalled with permanent fails
    my $setup   = shift;
    my $give_up = 3600;
    my $debug   = VRPipeTest::debug();
    $manager->set_verbose_global(1) if $debug;
    my $gave_up = 0;
    while (1) {
        last if setup_done($setup);
        
        if ($give_up-- <= 0) {
            $gave_up = 1;
            warn "setup ", $setup->id, " not finished yet, but giving up after 1500 cycles\n";
            last;
        }
        
        sleep(1);
    }
    
    $manager->set_verbose_global(0) if $debug;
    return $gave_up ? 0 : 1;
}

sub setup_done {
    my $setup = shift;
    
    # has it started?
    $setup->datasource->elements;
    my $found = VRPipe::DataElement->search({ datasource => $setup->datasource->id, withdrawn => 0 });
    return 0 unless $found;
    
    my $setup_id  = $setup->id;
    my $pipeline  = $setup->pipeline;
    my $num_steps = $pipeline->step_members;
    
    # have all dataelements completed all steps?
    my @not_complete = VRPipe::DataElementState->search({ pipelinesetup => $setup_id, 'dataelement.withdrawn' => 0, completed_steps => { '!=' => $num_steps } }, { join => 'dataelement' });
    return 1 if !@not_complete;
    
    # have all incomplete dataelements stalled due to failures?
    my $incomplete_due_to_stall = 0;
    foreach my $des (@not_complete) {
        my $de = $des->dataelement;
        my ($ss) = VRPipe::StepState->search({ complete => 0, dataelement => $de, pipelinesetup => $setup });
        next unless $ss;
        my $stalled_or_complete_subs = 0;
        next unless $ss->dispatched;
        my @subs = VRPipe::Submission->search({ stepstate => $ss });
        foreach my $sub (@subs) {
            if ($sub->_done) {
                $stalled_or_complete_subs++;
            }
            elsif ($sub->_failed && $sub->retries == 3) {
                $stalled_or_complete_subs++;
            }
        }
        if ($stalled_or_complete_subs == @subs) {
            $incomplete_due_to_stall++;
        }
    }
    
    return $incomplete_due_to_stall == @not_complete;
}

sub get_bam_header {
    my $bam      = shift;
    my $bam_path = $bam->absolute;
    my $header   = `samtools view -H $bam_path`;
    return split /\n/, $header;
}

sub get_bam_records {
    my $bam      = shift;
    my $bam_path = $bam->absolute;
    my $records  = `samtools view $bam_path`;
    return split /\n/, $records;
}

sub finish {
    exit;
}

1;
