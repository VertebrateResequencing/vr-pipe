
=head1 NAME

VRPipe::Manager - methods for managing the execution of pipelines

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This is the main module used to discover what work needs to be done ('trigger'
each L<VRPipe::PipelineSetup>) and then to 'dispatch' that work out to the
system's job scheduler to actually run the command lines on the compute
cluster.

The B<vrpipe-trigger_pipelines> and B<vrpipe-dispatch_pipelines> scripts call
methods of this module.

Note that the current implementation is slow and inefficient, with lots of
serial looping over thousands of objects. A radical overhaul is planned.

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

use VRPipe::Base;

class VRPipe::Manager extends VRPipe::Persistent {
    use Parallel::ForkManager;
    use Sys::CPU;
    use POSIX qw(ceil);
    
    our $DEFAULT_MAX_PROCESSES = Sys::CPU::cpu_count();
    our %step_limits;
    our %setups_with_step_limits;
    our %setups_checked_for_step_limits;
    our %step_limit_instigators;
    our $do_step_limits = 0;
    
    has '_running' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has '_run_start' => (
        is          => 'rw',
        isa         => Datetime,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has '_instance_run_start' => (
        is  => 'rw',
        isa => Datetime
    );
    
    has 'global_limit' => (
        is  => 'rw',
        isa => PositiveInt
    );
    
    # public getters for our private attributes
    method running {
        return $self->_running;
    }
    
    method run_start {
        return $self->_run_start;
    }
    
    __PACKAGE__->make_persistent();
    
    around get (ClassName|Object $self: Persistent :$id?, PositiveInt :$global_limit = 500) {
        return $self->create(global_limit => $global_limit);
    }
    
    around create (ClassName|Object $self: Persistent :$id?, PositiveInt :$global_limit = 500) {
        # our db-storage needs are class-wide, so we only have 1 row in our
        # table
        return $self->$orig(id => 1, global_limit => $global_limit);
    }
    
    method register_farm_server (Str $farm, Bool :$only_ours = 0) {
        my $transaction = sub {
            my ($fs) = VRPipe::FarmServer->search({ farm => $farm }, { for => 'update' });
            return if ($fs && $fs->alive);
            return VRPipe::FarmServer->create(farm => $farm, only_ours => $only_ours);
        };
        
        return $self->do_transaction($transaction, "Could not register farm $farm");
    }
    
    method setups (Str :$pipeline_name?) {
        my @setups;
        foreach my $ps (VRPipe::PipelineSetup->search({}, { prefetch => 'pipeline' })) {
            my $p = $ps->pipeline;
            if ($pipeline_name) {
                next unless $p->name eq $pipeline_name;
            }
            
            # we can have step-specific limits on how many subs to run at
            # once; see if any are in place and set them globally for all
            # setups
            if ($ps->active && !exists $setups_checked_for_step_limits{$ps_id}) {
                my $user_opts = $ps->options;
                
                my @stepms = $ps->pipeline->step_members;
                foreach my $stepm (@stepms) {
                    my $step      = $stepm->step;
                    my $step_name = $step->name;
                    
                    my $limit = $user_opts->{ $step_name . '_max_simultaneous' };
                    unless ($limit) {
                        $limit = $step->max_simultaneous;
                    }
                    
                    if ($limit) {
                        if (!defined $step_limits{$step_name} || $limit <= $step_limits{$step_name}) {
                            $step_limits{$step_name}                                = $limit;
                            $setups_with_step_limits{$ps_id}->{$step_name}          = $limit;
                            $step_limit_instigators{$step_name}->{$limit}->{$ps_id} = 1;
                            $do_step_limits                                         = 1;
                        }
                    }
                }
                
                $setups_checked_for_step_limits{$ps_id} = 1;
            }
            elsif (!$ps->active && exists $setups_with_step_limits{$ps_id}) {
                my @limited_steps = keys %{ $setups_with_step_limits{$ps_id} };
                while (my ($step_name, $limit) = each %{ $setups_with_step_limits{$ps_id} }) {
                    delete $step_limit_instigators{$step_name}->{$limit}->{$ps_id};
                    my $current_limit = $step_limits{$step_name};
                    unless (keys %{ $step_limit_instigators{$step_name}->{$current_limit} }) {
                        delete $step_limit_instigators{$step_name}->{$current_limit};
                        my @other_limits = sort { $a <=> $b } keys %{ $step_limit_instigators{$step_name} };
                        if (@other_limits) {
                            $step_limits{$step_name} = $other_limits[0];
                        }
                        else {
                            delete $step_limits{$step_name};
                            delete $step_limit_instigators{$step_name};
                        }
                    }
                }
                
                delete $setups_with_step_limits{$ps_id};
                delete $setups_checked_for_step_limits{$ps_id};
                
                unless (keys %setups_with_step_limits) {
                    $do_step_limits = 0;
                }
            }
            
            push(@setups, $ps);
        }
        
        return @setups;
    }
    
    method handle_submissions (Int :$max_retries = 3) {
        # call setups to find step limits
        $self->setups;
        
        my $unfinished = VRPipe::Submission->search({ '_done' => 0 }) || 0;
        $self->debug("There are $unfinished unfinished submissions to work on");
        
        if ($unfinished) {
            $self->check_running();
            $self->resubmit_failures(max_retries => $max_retries);
            $self->batch_and_submit();
            return 0;
        }
        else {
            return 1;
        }
    }
    
    method check_running {
        # update the status of each submission in case any of them finished
        my $pager = VRPipe::Submission->search_paged({ '_done' => 0, '_failed' => 0 }, { prefetch => 'job' });
        
        my $still_not_done = 0;
        my $c              = 0;
        while (my $subs = $pager->next) {
            foreach my $sub (@$subs) {
                my $job = $sub->job;
                $c++;
                $self->debug("loop $c, sub " . $sub->id . " job " . $job->id);
                if ($job->running) {
                    # check we've had a recent heartbeat
                    $self->debug(" -- running...");
                    if ($job->unresponsive) {
                        $self->debug(" -- unresponsive");
                        $job->kill_job;
                        $sub->update_status();
                    }
                    elsif ($sub->close_to_time_limit(30)) {
                        # user's scheduler might kill the submission if it runs too
                        # long in the queue it was initially submitted to; avoid
                        # this by changing queue as necessary
                        $self->debug(" -- within 30mins of time limit, will increase by 2hrs...");
                        $sub->extra_time(2);
                    }
                }
                elsif ($job->finished) {
                    $self->debug(" -- finishing by calling sub->update_status...");
                    $sub->update_status();
                    $self->debug(" -- success");
                    next if $sub->done;
                }
                elsif ($sub->scheduled) {
                    # we may be pending in the scheduler, but if we've been
                    # apparently pending for ages, check the scheduler really is
                    # pending us and if not reset the submission
                    my $pend_time = $sub->pend_time;
                    $self->debug(" -- pending for $pend_time seconds...");
                    if ($pend_time > 21600) { #*** rather than check every time check_running is called after 6hrs have passed, can we only check again once every subsequent hour?
                        $self->warn("submission " . $sub->id . " has been scheduled for $pend_time seconds - will try and resolve this");
                        $sub->unschedule_if_not_pending;
                    }
                }
                
                $still_not_done++;
            }
        }
        
        $self->debug("There are $still_not_done submissions that still aren't done");
    }
    
    method resubmit_failures (Int :$max_retries) {
        # this checks all failed submissions, and uses the standard Submission
        # and Job methods to work out why and resubmit them as appropriate,
        # potentially with updated Requirements. max_retries defaults to 3.
        my $pager = VRPipe::Submission->search_paged({ _failed => 1, retries => { '<' => $max_retries } }, { prefetch => [qw(scheduler requirements)] });
        
        while (my $subs = $pager->next) {
            foreach my $sub (@$subs) {
                my $parser = $sub->scheduler_stdout;
                if ($parser) {
                    my $status = $parser->status;
                    #*** how can we make parser status method be present for all
                    #    scheduler outputs, and how can we make the return
                    #    values generic?
                    if ($status) {
                        if ($status eq 'MEMLIMIT') {
                            $sub->extra_memory;
                        }
                        elsif ($status eq 'RUNLIMIT') {
                            my $hrs         = ceil($parser->time / 60 / 60) + 1;
                            my $current_hrs = $sub->time;
                            $sub->extra_time($hrs - $current_hrs);
                        }
                    }
                }
                $sub->retry;
            }
        }
    }
    
    method batch_and_submit {
        # this groups all the non-permanently-failed ->submissions into lists
        # based on shared Requirements objects, creates [PersistentArray,
        # Requirements] for each group, and uses those to do a Scheduler->submit
        # for each group. If a group consists of only 1 Submission, then it
        # won't bother with the creating a PersistentArray step.
        
        # if any step limits are in place we first have to see how many of each
        # step are currently in the scheduler
        my %step_counts;
        if ($do_step_limits) {
            my $array_ref = VRPipe::Submission->get_column_values(
                ['step.name', 'count(step.name)'],
                { '_done' => 0, '_failed' => 0, '_sid' => { '!=', undef } },
                {
                    join     => { stepstate => { stepmember => 'step' } },
                    group_by => ['step.name']
                }
            );
            foreach my $vals (@$array_ref) {
                $step_counts{ $vals->[0] } = $vals->[1];
            }
        }
        
        # we want to limit the number of jobs we have running based on the
        # global limit, so we'll count the number of submissions we batch and
        # end the submission loop when over the limit. We'll also start the
        # count at the number of currently running jobs.
        my $limit = $self->global_limit;
        my $count = VRPipe::Job->search({ 'running' => 1 });
        
        # we want to now find all submissions that are not done or failed and
        # that haven't already been submitted.
        my $pager = VRPipe::Submission->search_paged({ '_done' => 0, '_failed' => 0, '_sid' => undef }, { order_by => 'requirements', prefetch => [qw(job requirements)] }, 1000);
        
        my $scheduler = VRPipe::Scheduler->get;
        SLOOP: while (my $subs = $pager->next) {
            my %batches;
            foreach my $sub (@$subs) {
                # if the job is a block_and_skip_if_ok job, we don't actually
                # block because of race condition issues, and because of fail
                # and restart issues. Instead we always only actually submit if
                # this submission is the first incomplete submission created for
                # this $job. This also saves us creating lots of submissions for
                # a single job, which is confusing and wasteful.
                my $job = $sub->job;
                if ($job->block_and_skip_if_ok) {
                    my ($first_sub) = VRPipe::Submission->search({ 'job' => $job->id, '_done' => 0 }, { rows => 1, order_by => { -asc => 'id' } });
                    next unless $first_sub->id == $sub->id;
                }
                
                # step limit handling
                if ($do_step_limits) {
                    my $step_name = $sub->stepstate->stepmember->step->name;
                    if (exists $step_limits{$step_name}) {
                        $step_counts{$step_name}++;
                        next if $step_counts{$step_name} > $step_limits{$step_name};
                    }
                }
                
                # global limit handling
                $count++;
                last SLOOP if $count >= $limit;
                
                # we're good to batch this one
                push(@{ $batches{ $sub->requirements->id } }, $sub);
            }
            
            while (my ($req_id, $subs) = each %batches) {
                $self->debug("_~_ Batched and submitted an array of " . scalar(@$subs) . " subs");
                my $sid = $scheduler->submit(
                    array        => $subs,
                    requirements => $subs->[0]->requirements
                );
            }
        }
    }
}

1;
