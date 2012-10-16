
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
        my $job = $self->job;
        $self->reselect_values_from_db;
        if ($self->_claim) {
            if ($job->alive) {
                return 0;
            }
            else {
                if ($job->heartbeat) {
                    # another process probably claimed this sub, started running
                    # the job, but then died without releasing the claim: just
                    # go ahead and take the claim for ourselves
                    $self->_reset_job;
                    $self->_own_claim(1);
                    return 1;
                }
                else {
                    # however, the other process may have claimed it and not yet
                    # called run() and started the job's heartbeat
                    return 0;
                    
                    #*** but we have an edge-case hole where the sub got claimed
                    #    but that process never called run() or release(), in
                    #    which case, we will always return 0 here and the job
                    #    will never get run!
                }
            }
        }
        else {
            if ($job->alive) {
                $self->_reset_job;
            }
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
        }
        else {
            $self->_done(0);
            $self->_failed(1);
        }
        
        $self->archive_output;
        
        $self->update;
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
                $self->throw("toa $count had dest $dest compared to source " . $source->path . "; should have been $path for vrfile $vrfile_id");
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
    
    method extra_time (Int $extra = 3600) {
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
        $file->s || return;
        return VRPipe::Parser->create($type, { file => $file->path }); # in case $file is not type $type, we send the path, not the object
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
        if ($job->start_time && !$job->end_time) {
            $job->kill_job;
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
