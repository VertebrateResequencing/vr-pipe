
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
    use DateTime;
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
    
    has '_own_claim' => (
        is      => 'rw',
        isa     => 'Bool',
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
    method claim {
        if ($self->_claim) {
            return $self->_own_claim ? 1 : 0;
        }
        else {
            $self->_claim(1);
            $self->update;
            $self->_own_claim(1);
            return 1;
        }
    }
    
    method release {
        return unless $self->_own_claim;
        $self->_claim(0);
        $self->_own_claim(0);
        $self->update;
    }
    
    method submit {
        $self->scheduler->submit(submission => $self);
    }
    
    around scheduler {
        if ($self->scheduled) {
            return $self->$orig;
        }
        else {
            return $self->$orig(@_);
        }
    }
    
    method update_status {
        $self->reselect_values_from_db;
        $self->throw("Cannot call update_status when the job " . $self->job->id . " is not finished") unless $self->job->finished;
        return if $self->done || $self->failed;
        
        if ($self->job->ok) {
            $self->_done(1);
            $self->_failed(0);
        }
        else {
            $self->_done(0);
            $self->_failed(1);
        }
        
        unless ($self->_hid) {
            # we're a submission for a job that completed in a different
            # submission, so we didn't actually run and have no output
            $self->update;
            return;
        }
        
        #*** these 2 calls are probably the cause of massive delays when going
        #    from jobs being finished to submissions being done... can we
        #    optimise?
        #    It is actualy the archive_output that causes the ~1s delays,
        #    and it is not the filesystem operations, but the retrieval of
        #    the 8 VRPipe::File objects... at the very least path column of file
        #    table must be indexed, or delays go up to ~6s+
        my $t1 = time();
        $self->sync_scheduler;
        my $t2 = time();
        $self->archive_output;
        my $t3 = time();
        $self->_sid(undef);
        if ($t3 - $t1 > 2) {
            my $ss_t = $t2 - $t1;
            my $ao_t = $t3 - $t2;
            $self->debug("in update_status, sync_scheduler took $ss_t seconds; archive_output took $ao_t seconds");
        }
        
        $self->update;
    }
    
    method sync_scheduler {
        my $sid = $self->sid;
        unless ($sid) {
            return;
        }
        
        # assume that if its been more than a minute since the job ended, the
        # scheduler must have generated its output by now, so we don't have to
        # ask the scheduler about this sid
        my $job      = $self->job;
        my $end_time = $self->job->end_time;
        if ($end_time && (time() - $end_time->epoch > 60)) {
            return;
        }
        
        $self->scheduler->wait_for_sid($sid, $self->_aid, 5);
    }
    
    method archive_output {
        my @to_archive;
        push(@to_archive, [$self->job->stdout_file, $self->job_stdout_file, 1]);
        push(@to_archive, [$self->job->stderr_file, $self->job_stderr_file, 1]);
        push(@to_archive, [$self->scheduler_stdout_file(orig => 1), $self->scheduler_stdout_file, 0]);
        push(@to_archive, [$self->scheduler_stderr_file(orig => 1), $self->scheduler_stderr_file, 1]);
        
        my $count = 0;
        foreach my $toa (@to_archive) {
            $count++;
            my ($source, $dest, $add_marker) = @$toa;
            next unless ($source && $dest);
            if (!ref($dest)) {
                $self->throw("toa $count had dest $dest compared to source " . $source->path);
            }
            $self->concatenate(
                $source, $dest,
                unlink_source => 1,
                add_marker    => $add_marker,
                max_lines     => 1000
            );
        }
    }
    
    # requirement passthroughs and extra_* methods
    method _add_extra (Str $type, Int $extra) {
        my $new_req = $self->requirements->clone($type => $self->$type() + $extra);
        
        # we need to check with the scheduler if the new requirements would put
        # us in a different queue, and if so switch queues on any submission
        # that is currently running the job
        my ($switch_queue, $scheduler);
        my $job = $self->job;
        if ($job->running) {
            $scheduler = $self->scheduler;
            my $current_queue = $scheduler->determine_queue($self->requirements);
            my $new_queue     = $scheduler->determine_queue($new_req);
            if ($new_queue ne $current_queue) {
                $switch_queue = $new_queue;
            }
        }
        
        # we want to add extra * for all submissions that are for this sub's
        # job, incase it is not this sub in an array of block_and_skip jobs that
        # gets retried; also, if we're switching queues, it will be the first
        # sub that is actually running the job
        my $pager = VRPipe::Submission->search_paged({ 'job' => $job->id }, { order_by => { -asc => 'id' } });
        my $first = 1;
        while (my $to_extras = $pager->next) {
            foreach my $to_extra (@$to_extras) {
                if ($first && $switch_queue && $to_extra->sid && !$to_extra->done) {
                    $self->debug("calling switch_queues(" . $to_extra->sid . ", $switch_queue)");
                    $scheduler->switch_queue($to_extra->sid, $switch_queue);
                }
                $first = 0;
                
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
        $self->_add_extra('memory', $extra);
    }
    
    method time {
        return $self->requirements->time;
    }
    
    method extra_time (Int $extra = 1) {
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
    
    sub DEMOLISH {
        return if in_global_destruction;
        my $self = shift;
        $self->release if $self->in_storage;
    }
    
    __PACKAGE__->make_persistent();
    
    method close_to_time_limit (Int $minutes = 30) {
        my $seconds = $minutes * 60;
        my $job     = $self->job;
        return $job->wall_time > (($self->requirements->time * 60 * 60) - $seconds);
    }
    
    # where have our scheduler and job stdout/err files gone?
    method std_dir (Bool $orig = 0) {
        my ($for) = $self->_for($orig);
        $for || return;
        return $self->scheduler->output_dir($for);
    }
    
    method _for (Bool $orig = 0) {
        my $hid = $self->_hid || return;
        my $for = $self->job;
        my $index;
        if ($self->_aid) {
            $for = VRPipe::PersistentArray->get(id => $hid) if $orig;
            $index = $self->_aid;
        }
        else {
            $for = $self if $orig;
        }
        return ($for, $index);
    }
    
    method _scheduler_std_file (Str $method where {$_ eq 'scheduler_output_file' || $_ eq 'scheduler_error_file'}, Str $type, Bool $orig = 0) {
        my $std_dir = $self->std_dir($orig) || return;
        my $std_io_file;
        if ($orig) {
            $std_io_file = $self->scheduler->$method($std_dir);
            my (undef, $index) = $self->_for($orig);
            if ($index) {
                $std_io_file .= '.' . $index;
            }
        }
        else {
            $std_io_file = file($std_dir, $method);
        }
        
        return VRPipe::File->create(path => $std_io_file, type => $type);
    }
    
    method scheduler_stdout_file (Bool :$orig = 0) {
        return $self->_scheduler_std_file('scheduler_output_file', substr($self->scheduler->type, 0, 3), $orig);
    }
    
    method scheduler_stderr_file (Bool :$orig = 0) {
        return $self->_scheduler_std_file('scheduler_error_file', 'cat', $orig);
    }
    
    method _std_parser (VRPipe::File $file, Str $type) {
        $file->s || return;
        return VRPipe::Parser->create($type, { file => $file->path }); # in case $file is not type $type, we send the path, not the object
    }
    
    method scheduler_stdout {
        my $file = $self->scheduler_stdout_file || return;
        unless ($file->s) {
            $file = $self->scheduler_stdout_file(orig => 1);
        }
        return $self->_std_parser($file, lc(substr($self->scheduler->type, 0, 3)));
    }
    
    method scheduler_stderr {
        my $file = $self->scheduler_stderr_file || return;
        unless ($file->s) {
            $file = $self->scheduler_stderr_file(orig => 1);
        }
        return $self->_std_parser($file, 'cat');
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
    
    method pend_time {
        my $job = $self->job;
        my $scheduled_time = $self->_scheduled || $self->throw("called pend_time, yet submission " . $self->id . " has not been scheduled!");
        $scheduled_time = $scheduled_time->epoch;
        if ($job->running || $job->finished) {
            return $job->start_time->epoch - $scheduled_time;
        }
        else {
            return time() - $scheduled_time;
        }
    }
    
    method unschedule_if_not_pending {
        my $sid       = $self->sid || return;
        my $aid       = $self->_aid;
        my $scheduler = $self->scheduler;
        
        my $status = $scheduler->sid_status($sid, $aid);
        return if $status eq 'PEND';
        $scheduler->kill_sid($sid, $aid, 5);
        
        $self->_reset;
    }
    
    method retry {
        return unless ($self->done || $self->failed);
        
        # reset the job
        $self->_reset_job;
        
        # reset ourself
        my $retries = $self->retries;
        $self->retries($retries + 1);
        $self->_reset;
    }
    
    method start_over {
        # reset the job
        $self->_reset_job;
        
        # reset ourself and also set retries to 0
        $self->retries(0);
        $self->_reset;
    }
    
    method _reset_job {
        my $job = $self->job;
        unless ($job->finished) {
            $job->kill_job;
        }
        
        my $ofiles = $job->output_files;
        foreach my $file (@$ofiles) {
            $file->unlink;
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
