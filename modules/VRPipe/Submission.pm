
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

Copyright (c) 2011-2013 Genome Research Limited.

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
    
    method claim_and_run (PositiveInt :$allowed_time?, Object :$redis?) {
        # we'll return one of a number of responses: 0 = there was some problem
        # so the job and sub were reset; 1 = this process just now claimed the
        # sub and started running the job; 2 = the job is already exited and the
        # sub marked as done/failed; 3 = another process is currently running
        # the job
        my $response;
        
        # we'll also return the job instance that we called run() on, since that
        # will hold non-persistent data for _signalled_to_death
        my $job_instance;
        my $transaction = sub {
            # lock the Sub and Job rows before trying to claim
            my $job = $self->job;
            $self->lock_row($job);
            $self->lock_row($self);
            
            my $ss             = $self->stepstate;
            my $ps             = $ss->pipelinesetup;
            my %log_event_args = (dataelement => $ss->dataelement->id, stepstate => $ss->id, submission => $self->id, job => $job->id);
            
            # first check if the job has already started or finished
            if ($self->done || $self->failed || (defined $job->exit_code && $job->end_time)) {
                if (!($self->done || $self->failed)) {
                    # the job ran under a different submission; we just need to
                    # update this submission
                    if ($job->ok) {
                        $self->_done(1);
                        $self->_failed(0);
                    }
                    else {
                        $self->_done(0);
                        $self->_failed(1);
                    }
                    $self->_claim(0);
                    $self->update;
                }
                elsif ($self->failed) {
                    $self->retry if $self->retries < 3;
                }
                $response = 2;
                $ps->log_event("claim_and_run() found that the Submission was already done or failed", %log_event_args);
                return;
            }
            elsif ($job->start_time) {
                if ($job->alive(no_suicide => 1)) {
                    $response = 3;
                    $ps->log_event("claim_and_run() found that the Job had already started and is currently alive", %log_event_args);
                    return;
                }
                else {
                    # looks like the process that previously started running
                    # the job exited before updating the job's state; reset
                    $ps->log_event("claim_and_run() found that the Job had already started yet was dead, so it will be reset", %log_event_args);
                    $self->_reset_job;
                    $self->_reset;
                    $response = 0;
                    return;
                }
            }
            elsif ($self->_claim) {
                # this should be impossible, since sub->_claim and
                # job->start_time are only ever set together in the same
                # transaction
                $ps->log_event("claim_and_run() found that the Job had not started yet the Submission was claimed, which should be impossible; they will both be reset", %log_event_args);
                $self->_reset_job;
                $self->_reset;
                $response = 0;
                return;
            }
            
            # see if we can get the job to start running, and set claim if so
            my $run_response = $job->run(submission => $self, $allowed_time ? (allowed_time => $allowed_time) : (), $redis ? (redis => $redis) : ());
            if (!$run_response) {
                $ps->log_event("claim_and_run() tried to run() the Job but it was already running", %log_event_args);
                $response = 3;
            }
            elsif ($run_response == -1) {
                if ($self->done || $self->failed) {
                    $ps->log_event("claim_and_run() tried to run() the Job but it was already complete", %log_event_args);
                    $response = 2;
                }
                else {
                    $ps->log_event("claim_and_run() Job->run() call claimed the Job was complete, yet the Submission is not marked as done or failed, which should be impossible; both will be reset", %log_event_args);
                    $self->_reset_job;
                    $self->_reset;
                    $response = 0;
                }
            }
            elsif ($run_response == 1) {
                $self->_claim(1);
                $self->update;
                $response     = 1;
                $job_instance = $job;
            }
            
            my $claim = $self->_claim    || 0;
            my $st    = $job->start_time || 'n/a';
        };
        $self->do_transaction($transaction, "Failed when trying to claim and run");
        
        return ($response, $job_instance);
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
    
    method archive_output {
        $self->throw("Cannot call archive_output when the job " . $self->job->id . " is not finished") unless $self->job->end_time;
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
        my $file = VRPipe::File->create(path => file($std_dir, 'job_std' . $kind), type => 'cat');
        $file->update_stats_from_disc;
        return $file;
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
        $job->reselect_values_from_db;
        if ($job->start_time && !$job->end_time) {
            $self->stepstate->pipelinesetup->log_event("Submission->_reset_job() call will kill the currently running Job", dataelement => $self->stepstate->dataelement->id, stepstate => $self->stepstate->id, submission => $self->id, job => $self->job->id);
            $job->kill_job($self);
        }
        else {
            # delete its output
            my $ofiles = $job->output_files;
            if ($ofiles) {
                foreach my $file (@$ofiles) {
                    $file->unlink;
                }
            }
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
