
=head1 NAME

VRPipe::Submission - holds state for Jobs that need to be run

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A L<VRPipe::Job> holds a command line that needs to be executed, but to execute
it it must be passed to the system's job scheduler, which usually would like to
know the memory and time (etc.) requirements (L<VRPipe::Requirements>).

Submission associates Job with Requirements, and tracks scheduler-related
state. We have here the methods to retry a failed Job, and increase memory and
time requirements (associating with the Submission a new Requirements object)
if necessary.

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

class VRPipe::Submission extends VRPipe::Persistent {
    use Devel::GlobalDestruction;
    use VRPipe::Parser;
    use POSIX qw(ceil);
    
    has 'job' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::Job'
    );
    
    has 'stepstate' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::StepState'
    );
    
    has 'requirements' => (
        is       => 'rw',
        isa      => Persistent,
        coerce   => 1,
        required => 1,                                 # even though we're not a key
        traits   => ['VRPipe::Persistent::Attributes'],
        # handles => [qw(memory time cpu tmp_space local_space custom)]), *** doesn't work for some reason, and we need them read-only anyway
        belongs_to => 'VRPipe::Requirements'
    );
    
    has 'scheduler' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        required   => 1,
        builder    => '_build_default_scheduler',
        traits     => ['VRPipe::Persistent::Attributes'],
        belongs_to => 'VRPipe::Scheduler'
    );
    
    has 'retries' => (
        is      => 'rw',
        isa     => IntSQL [4],
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has '_done' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has '_failed' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has '_claim' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    method _build_default_scheduler {
        return VRPipe::Scheduler->create();
    }
    
    # public getters for our private attributes
    method done {
        return $self->_done;
    }
    
    method failed {
        return $self->_failed;
    }
    
    # other methods
    method _get_claim {
        my $claimed_by_us;
        my $transaction = sub {
            my ($sub_to_claim) = VRPipe::Submission->search({ id => $self->id }, { for => 'update' });
            my $job = $sub_to_claim->job;
            if ($sub_to_claim->_claim) {
                unless ($job->alive(no_suicide => 1)) {
                    # clean up any 'stuck' submissions
                    if ($job->heartbeat) {
                        # another process probably claimed this sub, started
                        # running the job, but then died without releasing the
                        # claim: clean this up now
                        $sub_to_claim->_claim(0);
                        $sub_to_claim->update;
                        $self->stepstate->pipelinesetup->log_event("claim_and_run() found that the Submission was claimed yet the Job was dead and had a heartbeat: the Job must have died so the claim was unset and will reset the Job", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $job->id);
                        $sub_to_claim->_reset_job;
                    }
                    else {
                        # however, the other process may have claimed it and not
                        # yet called run() and started the job's heartbeat
                        
                        # but we have an edge-case hole where the sub got
                        # claimed but that process never called run() or
                        # release(), in which case, we will always return 0 here
                        # and the job will never get run! So we wait 60s for
                        # the job to start running, and if it doesn't we release
                        # the claim
                        my $started_running = 0;
                        for (1 .. 12) {
                            sleep(5);
                            $job->reselect_values_from_db;
                            if ($job->alive(no_suicide => 1)) {
                                $started_running = 1;
                                last;
                            }
                            
                            #*** there's also some kind of issue that we can
                            # have no heartbeat yet ran the Job successfully
                            # already, and we might also have set Submission to
                            # done... not sure what to do about failures to
                            # update database as expected...
                        }
                        
                        unless ($started_running) {
                            $sub_to_claim->_claim(0);
                            $sub_to_claim->update;
                            $self->stepstate->pipelinesetup->log_event("claim_and_run() found that the Submission was claimed yet the Job was dead and had no heartbeat; after waiting 60s the Job still hasn't started so the claim was unset", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $job->id);
                            
                            unless (defined $job->exit_code && $job->end_time) {
                                $self->stepstate->pipelinesetup->log_event("claim_and_run() found that the Submission was claimed yet the Job was dead and had no heartbeat; after waiting 60s the Job had no end_time or exit_code so the Job will be reset", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $job->id);
                                $sub_to_claim->_reset_job;
                            }
                        }
                    }
                }
                
                # regardless of if it was stuck, we don't claim. If we fixed
                # being stuck, the next attempt to claim will work
                $claimed_by_us = 0;
            }
            else {
                # ok, looks like we're good to get the claim, but there's an
                # edge case where another process started running the job yet
                # database weirdness prevented the claim from being set to 1,
                # and then the process died before job could be updated; check
                # that the job isn't dead
                if ($job->alive(no_suicide => 1)) {
                    $sub_to_claim->_claim(1);
                    $sub_to_claim->update;
                }
                elsif ($job->start_time) {
                    $self->stepstate->pipelinesetup->log_event("claim_and_run() found that the Submission was unclaimed yet the Job is dead having started: will reset the Job", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $job->id);
                    $sub_to_claim->_reset_job;
                }
                else {
                    $sub_to_claim->_claim(1);
                    $sub_to_claim->update;
                    $claimed_by_us = 1;
                }
            }
        };
        $self->do_transaction($transaction, "Failed when trying to claim submission");
        
        $self->reselect_values_from_db;
        return $claimed_by_us;
    }
    
    method claim_and_run (PositiveInt :$allowed_time?) {
        my $run_ok;
        my $transaction = sub {
            if ($self->_get_claim) {
                if ($self->failed) {
                    $self->stepstate->pipelinesetup->log_event("claim_and_run() failed to get the claim, and since the Submission is failed will retry()", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
                    $self->retry;
                    $run_ok = 0;
                }
                else {
                    # start running the associated job
                    $self->stepstate->pipelinesetup->log_event("claim_and_run() got the claim so will call Job->run()", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
                    $run_ok = $self->job->run(submission => $self, $allowed_time ? (allowed_time => $allowed_time) : ());
                }
                
                unless ($run_ok && $self->job->start_time) {
                    $self->release;
                    $self->stepstate->pipelinesetup->log_event("claim_and_run() got the claim but Job->run did not result in the Job starting, so the claim was released", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
                }
            }
        };
        $self->do_transaction($transaction, "Failed when trying to claim and run");
        $self->disconnect;
        
        return $run_ok;
    }
    
    method release {
        $self->_claim(0);
        $self->update;
    }
    
    around scheduler {
        if ($self->_claim) {
            return $self->$orig;
        }
        else {
            return $self->$orig(@_);
        }
    }
    
    method update_status {
        $self->reselect_values_from_db;
        $self->throw("Cannot call update_status when the job " . $self->job->id . " is not finished") unless $self->job->end_time;
        return if $self->done || $self->failed;
        
        if ($self->job->ok) {
            $self->_done(1);
            $self->_failed(0);
            $self->stepstate->pipelinesetup->log_event("update_status() call found that the Job was ok, so done set to 1", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
        }
        else {
            $self->_done(0);
            $self->_failed(1);
            $self->stepstate->pipelinesetup->log_event("update_status() call found that the Job was not ok, so failed set to 1", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
        }
        $self->update;
        
        $self->archive_output;
    }
    
    method archive_output {
        my @to_archive;
        push(@to_archive, [$self->job->stdout_file, $self->job_stdout_file]);
        push(@to_archive, [$self->job->stderr_file, $self->job_stderr_file]);
        
        my $count = 0;
        foreach my $toa (@to_archive) {
            $count++;
            my ($source, $dest) = @$toa;
            next unless ($source && $dest);
            if (!ref($dest)) {
                my $std_dir = $self->std_dir || 'no_std_dir';
                my $path = file($std_dir, 'job_stderr');
                my $vrfile = VRPipe::File->create(path => $path, type => 'cat');
                my $vrfile_id = $vrfile ? $vrfile->id : 'no_id';
                $self->warn("toa $count had dest $dest compared to source " . $source->path . "; should have been $path for vrfile $vrfile_id");
                return;
            }
            $self->concatenate(
                $source, $dest,
                unlink_source => 1,
                add_marker    => 1,
                max_lines     => 1000
            );
        }
    }
    
    # requirement passthroughs and extra_* methods
    method _add_extra (Str $type, Int $extra) {
        my $new_req = $self->requirements->clone($type => $self->$type() + $extra);
        
        # we want to add extra * for all submissions that are for this sub's
        # job, incase it is not this sub in an array of block_and_skip jobs that
        # gets retried *** is this still possible?
        my $pager = VRPipe::Submission->search_paged({ 'job' => $self->job->id }, { order_by => { -asc => 'id' } });
        while (my $to_extras = $pager->next) {
            foreach my $to_extra (@$to_extras) {
                $to_extra->requirements($new_req);
                $to_extra->update;
            }
        }
        
        $self->reselect_values_from_db;
        return 1;
    }
    
    method memory {
        return $self->requirements->memory;
    }
    
    method extra_memory (Int $extra?) {
        unless ($extra) {
            # increase by 1GB or 30%, whichever is greater, and round up to
            # nearest 100
            #*** increase to greater than max seen in stepstats?
            my $minimum_memory_increase    = 1000;
            my $memory_increase_percentage = 0.3;
            my $current_mem                = $self->memory;
            my $updated_memory_limit       = $current_mem * (1 + $memory_increase_percentage);
            if ($updated_memory_limit < $current_mem + $minimum_memory_increase) {
                $updated_memory_limit = $current_mem + $minimum_memory_increase;
            }
            $extra = $updated_memory_limit - $current_mem;
        }
        
        $extra = ceil($extra / 100) * 100;
        
        $self->stepstate->pipelinesetup->log_event("extra_memory() call will change our Requirements by adding ${extra}MB", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
        $self->_add_extra('memory', $extra);
    }
    
    method time {
        return $self->requirements->time;
    }
    
    method extra_time (Int $extra = 3600) {
        $self->stepstate->pipelinesetup->log_event("extra_time() call will change our Requirements by adding ${extra}secs", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
        $self->_add_extra('time', $extra);
    }
    
    method cpu {
        return $self->requirements->cpu;
    }
    
    method extra_cpu (Int $extra = 1) {
        $self->_add_extra('cpu', $extra);
    }
    
    method tmp_space {
        return $self->requirements->tmp_space;
    }
    
    method extra_tmp_space (Int $extra = 4000) {
        $self->_add_extra('tmp_space', $extra);
    }
    
    method local_space {
        return $self->requirements->local_space;
    }
    
    method extra_local_space (Int $extra = 4000) {
        $self->_add_extra('local_space', $extra);
    }
    
    method custom {
        return $self->requirements->custom;
    }
    
    __PACKAGE__->make_persistent();
    
    method close_to_time_limit (Int $minutes = 5) {
        my $seconds = $minutes * 60;
        my $job     = $self->job;
        return $job->wall_time > (($self->requirements->time) - $seconds);
    }
    
    # where have our job stdout/err files gone?
    method std_dir (Bool $orig = 0) {
        my $for = $orig ? $self : $self->job;
        $for || return;
        return $self->scheduler->output_dir($for);
    }
    
    method _job_std_file (Str $kind where {$_ eq 'out' || $_ eq 'err'}) {
        # this is where we want the job stdo/e to be archived to, not where the
        # job initially spits it out to
        my $std_dir = $self->std_dir || return;
        return VRPipe::File->create(path => file($std_dir, 'job_std' . $kind), type => 'cat');
    }
    
    method job_stdout_file {
        return $self->_job_std_file('out');
    }
    
    method job_stderr_file {
        return $self->_job_std_file('err');
    }
    
    method job_stdout {
        my $file = $self->job_stdout_file || return;
        unless ($file->s) {
            $file = $self->job->stdout_file;
        }
        return $self->_std_parser($file, 'cat');
    }
    
    method job_stderr {
        my $file = $self->job_stderr_file || return;
        unless ($file->s) {
            $file = $self->job->stderr_file;
        }
        return $self->_std_parser($file, 'cat');
    }
    
    method _std_parser (VRPipe::File $file, Str $type) {
        # these files have a habit of disappearing?? Just check with disc that
        # they really exist
        $file->update_stats_from_disc;
        $file->s || return;
        return VRPipe::Parser->create($type, { file => $file->path }); # in case $file is not type $type, we send the path, not the object
    }
    
    method retry {
        return unless ($self->done || $self->failed);
        
        $self->stepstate->pipelinesetup->log_event("Submission->retry() call will reset our Job and the Submission", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
        
        # reset the job
        $self->_reset_job;
        
        # reset ourself
        my $retries = $self->retries;
        $self->retries($retries + 1);
        $self->_reset;
    }
    
    method start_over {
        my $ss = $self->stepstate;
        $ss->pipelinesetup->log_event("Submission->start_over was called", stepstate => $ss->id, dataelement => $ss->dataelement->id, submission => $self->id, job => $self->job->id, record_stack => 1);
        
        # reset the job
        $self->_reset_job;
        
        # reset ourself and also set retries to 0
        $self->retries(0);
        $self->_reset;
        
        $ss->pipelinesetup->log_event("Submission->start_over call returning", stepstate => $ss->id, dataelement => $ss->dataelement->id, submission => $self->id, job => $self->job->id);
    }
    
    method _reset_job {
        my $job = $self->job;
        if ($job->start_time && !$job->end_time) {
            $self->stepstate->pipelinesetup->log_event("Submission->_reset_job() call will kill the currently running Job", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
            $job->kill_job($self);
        }
        
        $job->reset_job;
    }
    
    method _reset {
        $self->_failed(0);
        $self->_done(0);
        $self->_claim(0);
        $self->update;
        $self->reselect_values_from_db;
    }
}

1;
