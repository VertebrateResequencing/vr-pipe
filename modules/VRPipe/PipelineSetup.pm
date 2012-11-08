
=head1 NAME

VRPipe::PipelineSetup - set up a pipeline

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The main thing that users want to do is "run a pipeline" on a given set of
data. Pipelines are modelled by L<VRPipe::Pipline> objects, and the set of data
is modelled by a L<VRPipe::DataSource>. A PipelineSetup relates the two and
also stores the user configuration of both, resulting in a named object
defining what the user wants to happen.

So users "set up" pipelines and get back a PipelineSetup. They, and the
B<VRPipe> system itself, look to these PipelineSetups to see what is supposed
to be run and how, on what data.

Multiple PipelineSetups can run at once, even working on the same Pipelines and
DataSources (presumably with different configuration options), and the system
ensures that there are no problems with similar work being done by different
PipelineSetups overwriting each other's files.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

use VRPipe::Base;

class VRPipe::PipelineSetup extends VRPipe::Persistent {
    use POSIX qw(ceil);
    
    has 'name' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'datasource' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::DataSource'
    );
    
    has 'pipeline' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::Pipeline'
    );
    
    has 'output_root' => (
        is     => 'rw',
        isa    => Dir,
        coerce => 1,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'options' => (
        is                   => 'rw',
        isa                  => 'HashRef',
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => sub { {} },
        allow_key_to_default => 1,
        is_key               => 1
    );
    
    has 'description' => (
        is          => 'rw',
        isa         => Text,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'active' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 1
    );
    
    has 'user' => (
        is      => 'rw',
        isa     => Varchar [64],
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 'vrpipe'
    );
    
    has 'desired_farm' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'controlling_farm' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    __PACKAGE__->make_persistent(has_many => [states => 'VRPipe::StepState']);
    
    # because lots of frontends get the dataelementstates for a particular
    # pipelinesetup
    sub _des_search_args {
        my $self              = shift;
        my $include_withdrawn = shift;
        my $only_withdrawn    = shift;
        my @terms             = ();
        if ($only_withdrawn) {
            @terms = ('dataelement.withdrawn' => 1);
        }
        elsif (!$include_withdrawn) {
            @terms = ('dataelement.withdrawn' => 0);
        }
        return ({ pipelinesetup => $self->id, @terms }, { prefetch => 'dataelement' });
    }
    
    method dataelementstates_pager (Bool :$include_withdrawn = 0, Bool :$only_withdrawn = 0) {
        return VRPipe::DataElementState->search_paged($self->_des_search_args($include_withdrawn, $only_withdrawn));
    }
    
    method dataelementstates (Bool :$include_withdrawn = 0, Bool :$only_withdrawn = 0) {
        return VRPipe::DataElementState->search($self->_des_search_args($include_withdrawn, $only_withdrawn));
    }
    
    around desired_farm (Maybe[Str] $farm?) {
        my $current_farm = $self->$orig();
        if ($farm) {
            if ($current_farm && $farm ne $current_farm) {
                $self->controlling_farm(undef);
            }
            return $self->$orig($farm);
        }
        return $current_farm;
    }
    
    method currently_complete {
        # we consider ourselves incomplete if we have no dataelementstates
        # assigned to us
        my $num_elements = VRPipe::DataElementState->search({ pipelinesetup => $self->id });
        return 0 unless $num_elements;
        
        # and we're incomplete if any of them have not completed all steps in
        # our pipeline
        my $num_steps_to_complete = $self->pipeline->step_members;
        my $elements_incomplete = VRPipe::DataElementState->search({ pipelinesetup => $self->id, completed_steps => { '<' => $num_steps_to_complete }, 'dataelement.withdrawn' => 0 }, { join => 'dataelement' });
        return $elements_incomplete ? 0 : 1;
    }
    
    method trigger (VRPipe::DataElement :$dataelement?, VRPipe::Interface::BackEnd :$backend?) {
        my $setup_id     = $self->id;
        my $pipeline     = $self->pipeline;
        my @step_members = $pipeline->step_members;
        my $num_steps    = scalar(@step_members);
        
        my $datasource  = $self->datasource;
        my $output_root = $self->output_root;
        $self->make_path($output_root);
        
        # we either loop through all incomplete elementstates, or the
        # (single) elementstate for the supplied dataelement
        my $pager;
        if ($dataelement) {
            $pager = VRPipe::DataElementState->search_paged({ pipelinesetup => $setup_id, dataelement => $dataelement->id, completed_steps => { '<', $num_steps }, 'dataelement.withdrawn' => 0 }, { prefetch => 'dataelement' });
        }
        else {
            $pager = $datasource->incomplete_element_states($self, prepare => 1);
        }
        my $all_done = 1;
        while (my $estates = $pager->next) {
            foreach my $estate (@$estates) {
                my $element         = $estate->dataelement;
                my $completed_steps = $estate->completed_steps;
                next if $completed_steps == $num_steps;
                
                my %previous_step_outputs;
                foreach my $member (@step_members) {
                    my $step_number = $member->step_number;
                    my $state       = VRPipe::StepState->create(
                        stepmember    => $member,
                        dataelement   => $element,
                        pipelinesetup => $self
                    );
                    
                    my $step = $member->step(previous_step_outputs => \%previous_step_outputs, step_state => $state);
                    if ($state->complete) {
                        $self->_complete_state($step, $state, $step_number, $pipeline, \%previous_step_outputs, $estate);
                        next;
                    }
                    
                    # have we previously done the dispatch dance and are
                    # currently waiting on submissions to complete?
                    my @submissions = $state->submissions;
                    if (@submissions) {
                        my $unfinished = VRPipe::Submission->search({ '_done' => 0, stepstate => $state->submission_search_id });
                        unless ($unfinished) {
                            my $ok = $step->post_process();
                            if ($ok) {
                                # we just completed all the submissions from a previous parse
                                $self->_complete_state($step, $state, $step_number, $pipeline, \%previous_step_outputs, $estate);
                                next;
                            }
                            else {
                                # we warn instead of throw, because the step may
                                # have discovered its output files are missing
                                # and restarted itself
                                $self->warn("submissions completed, but post_process failed");
                            }
                        }
                        # else we have $unfinished unfinished submissions from a
                        # previous parse and are still running
                    }
                    else {
                        # this is the first time we're looking at this step for
                        # this data member for this pipelinesetup
                        my $completed;
                        try {
                            $completed = $step->parse();
                        }
                        catch ($err) {
                            if ($backend) {
                                my $mt = VRPipe::MessageTracker->create(subject => "overall state of setup $setup_id");
                                unless ($mt->already_sent("parsing problem")) {
                                    $backend->log("When trying to parse step " . $step->name . " for setup $setup_id we hit the following error:\n$err", email_to => [$self->user], email_admin => 1, subject => "Setup $setup_id has problems");
                                }
                            }
                            else {
                                $self->verbose(1);
                                $self->warn("There is a problem with setup $setup_id (and we had no backend to send an email about this):\n" . $err);
                                $self->verbose(0);
                            }
                            $all_done = 0;
                            last;
                        }
                        
                        if ($completed) {
                            # on instant complete, parse calls post_process
                            # itself and only returns true if that was
                            # successfull
                            $self->_complete_state($step, $state, $step_number, $pipeline, \%previous_step_outputs, $estate);
                            next;
                        }
                        else {
                            my $dispatched = $step->dispatched();
                            if (@$dispatched) {
                                # is there another stepstate that we already
                                # made equivalent submissions for?
                                my %other_states;
                                foreach my $arrayref (@$dispatched) {
                                    my ($cmd, undef, $job_args) = @$arrayref;
                                    my ($job) = VRPipe::Job->search({ dir => $output_root, $job_args ? (%{$job_args}) : (), cmd => $cmd });
                                    last unless $job;
                                    my @submissions = VRPipe::Submission->search({ job => $job->id }, { prefetch => 'stepstate' });
                                    last unless @submissions;
                                    foreach my $sub (@submissions) {
                                        $other_states{ $sub->stepstate->id }++;
                                    }
                                }
                                delete $other_states{ $state->id };
                                
                                my $same_as_us;
                                my $needed_count = scalar @$dispatched;
                                foreach my $other_id (sort { $a <=> $b } keys %other_states) {
                                    my $count = $other_states{$other_id};
                                    next unless $needed_count == $count;
                                    $same_as_us = $other_id;
                                    last;
                                }
                                
                                if ($same_as_us) {
                                    # we just say that $state's submissions are
                                    # the same as the other other stepstate's
                                    $state->same_submissions_as($same_as_us);
                                    $state->update;
                                }
                                else {
                                    # create new submissions for the relevant
                                    # stepstate ($state may have had start_over
                                    # run on it, which would have deleted its
                                    # same_submissions_as stepstate's subs, and
                                    # we want to create the new submissions for
                                    # that same_submissions_as stepstate, not
                                    # for $state, hence the use of
                                    # $state->submission_search_id)
                                    foreach my $arrayref (@$dispatched) {
                                        my ($cmd, $reqs, $job_args) = @$arrayref;
                                        my $sub = VRPipe::Submission->create(job => VRPipe::Job->create(dir => $output_root, $job_args ? (%{$job_args}) : (), cmd => $cmd), stepstate => $state->submission_search_id, requirements => $reqs);
                                    }
                                }
                            }
                            else {
                                # it is possible for a parse to result in a
                                # different step being started over because
                                # input files were missing
                                $self->debug("step " . $step->id . " for data element " . $element->id . " for pipeline setup " . $self->id . " neither completed nor dispatched anything!");
                            }
                        }
                    }
                    
                    $all_done = 0;
                    last;
                }
            }
        }
        
        return $all_done;
    }
    
    method _complete_state (VRPipe::Step $step, VRPipe::StepState $state, Int $step_number, VRPipe::Pipeline $pipeline, PreviousStepOutput $previous_step_outputs, VRPipe::DataElementState $estate) {
        while (my ($key, $val) = each %{ $step->outputs() }) {
            $previous_step_outputs->{$key}->{$step_number} = $val;
        }
        unless ($state->complete) {
            unless ($state->same_submissions_as) {
                # are there any behaviours to trigger?
                foreach my $behaviour (VRPipe::StepBehaviour->search({ pipeline => $pipeline->id, after_step => $step_number })) {
                    $behaviour->behave(data_element => $state->dataelement, pipeline_setup => $state->pipelinesetup);
                }
                
                # add to the StepStats
                my %done_jobs;
                foreach my $submission ($state->submissions) {
                    my $job = $submission->job;
                    next if $done_jobs{ $job->id };
                    VRPipe::StepStats->create(step => $step, pipelinesetup => $state->pipelinesetup, submission => $submission, memory => $job->peak_memory || 0, time => $job->wall_time);
                    $done_jobs{ $job->id } = 1;
                }
            }
            
            $state->complete(1);
            $state->update;
        }
        
        my $completed_steps = $estate->completed_steps;
        if ($step_number > $completed_steps) {
            $estate->completed_steps($step_number);
            $estate->update;
        }
    
    }
}

1;
