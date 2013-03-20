
=head1 NAME

VRPipe::LocalScheduler - a job scheduler for the local node

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

For testing purposes, or for those without a compute cluster, this simple job
scheduler can be used instead of LSF et al. It implements a single queue of
jobs, and it just tries to blindly execute them on the local CPU, without
regard for their memory needs or time requirements.

Note that it is poorly tested and perhaps unreliable, so probably should not be
used in production or at scale.

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

class VRPipe::LocalScheduler {
    use AnyEvent::ForkManager;
    use Sys::CPU;
    use Cwd;
    use POSIX qw(ceil);
    
    our $DEFAULT_CPUS = Sys::CPU::cpu_count();
    
    has 'cpus' => (
        is      => 'ro',
        isa     => 'Int',
        default => $DEFAULT_CPUS
    );
    
    # for submit
    has 'o' => (
        is            => 'ro',
        isa           => 'Str',
        documentation => 'absolute path to stdout of scheduler (required for submit)'
    );
    
    has 'e' => (
        is            => 'ro',
        isa           => 'Str',
        documentation => 'absolute path to stderr of scheduler (required for submit)'
    );
    
    has 'a' => (
        is            => 'ro',
        isa           => 'Int',
        documentation => 'to specify a job array give the size of the array (default 1)',
        default       => 1
    );
    
    method submit (Str $cmd, HashRef $env) {
        my $o_file = $self->o;
        my $e_file = $self->e;
        unless ($o_file && $e_file) {
            $self->throw("-o and -e are required for submit");
        }
        
        my $array_size = $self->a;
        my $lsj = VRPipe::LocalSchedulerJob->create(cmd => $cmd, array_size => $array_size, cwd => cwd(), env => $env);
        
        my $user = getlogin || getpwuid($<);
        my @lsjs_args;
        foreach my $aid (1 .. $array_size) {
            my $this_o = $o_file;
            $this_o =~ s/\%I/$aid/g;
            my $this_e = $e_file;
            $this_e =~ s/\%I/$aid/g;
            push(@lsjs_args, { localschedulerjob => $lsj, aid => $aid, o_file => $this_o, e_file => $this_e, user => $user });
        }
        VRPipe::LocalSchedulerJobState->bulk_create_or_update(@lsjs_args);
        
        my $sid = $lsj->id;
        return "Job <$sid> is submitted";
    }
    
    method jobs (ArrayRef $ids) {
        my @lsjss;
        my @warnings;
        if (@$ids) {
            foreach my $id (@$ids) {
                my ($job_states, $warning) = $self->job_states($id);
                if ($warning) {
                    push @warnings, $warning;
                }
                else {
                    push(@lsjss, @$job_states);
                }
            }
        }
        else {
            foreach my $lsj (VRPipe::LocalSchedulerJob->search({})) {
                push(@lsjss, $lsj->jobstates);
            }
        }
        
        unless (@lsjss) {
            return ([], \@warnings);
        }
        
        my @lines;
        push(@lines, join("\t", qw(JOBID INDEX STAT USER EXEC_HOST SUBMIT_TIME RUN_TIME CMD)));
        foreach my $lsjs (@lsjss) {
            my $lsj = $lsjs->localschedulerjob;
            push(@lines, join("\t", $lsj->id, $lsjs->aid, $lsjs->current_status, $lsjs->user, $lsjs->host || '', $lsj->creation_time, $lsjs->run_time, $lsj->cmd));
        }
        
        return (\@lines, \@warnings);
    }
    
    method job_states ($id) {
        my ($sid, $aid);
        if ($id =~ /(\d+)\[(\d+)\]/) {
            ($sid, $aid) = ($1, $2);
        }
        elsif ($id =~ /^\d+$/) {
            $sid = $id;
        }
        else {
            die "bad id format '$id'";
        }
        
        my ($lsj) = VRPipe::LocalSchedulerJob->search({ id => $sid });
        
        unless ($lsj) {
            return ([], "Job <$sid> is not found");
        }
        
        my @lsjss;
        foreach my $lsjs (VRPipe::LocalSchedulerJobState->search({ localschedulerjob => $sid })) {
            if ($aid) {
                next unless $lsjs->aid == $aid;
            }
            push(@lsjss, $lsjs);
        }
        if ($aid && !@lsjss) {
            return ([], "Job <$sid\[$aid\]> is not found");
        }
        
        return (\@lsjss);
    }
    
    method unfinished_lsjss {
        return VRPipe::LocalSchedulerJobState->search({ start_time => undef });
    }
    
    method process_queue {
        my @lsjss = $self->unfinished_lsjss;
        @lsjss || return;
        
        my $fm = AnyEvent::ForkManager->new(max_workers => $self->cpus);
        
        foreach my $lsjs (@lsjss) {
            $fm->start(
                cb => sub {
                    my ($fm, $lsjs) = @_;
                    $lsjs->start_job;
                },
                args => [$lsjs]
            );
        }
        
        $fm->wait_all_children;
        
        return 1;
    }
    
    method kill (ArrayRef $ids) {
        my @lsjss;
        my @warnings;
        my @lines;
        foreach my $id (@$ids) {
            my ($job_states, $warning) = $self->job_states($id);
            if ($warning) {
                push @warnings, $warning;
            }
            else {
                push(@lsjss, @$job_states);
                push @lines, "Job $id will be killed";
            }
        }
        
        foreach my $lsjss (@lsjss) {
            $lsjss->kill_job;
        }
        return (\@lines, \@warnings);
    }
}
