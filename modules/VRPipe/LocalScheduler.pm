=head1 NAME

VRPipe::LocalScheduler - a job scheduler for the local node

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

For testing purposes, or for those without a compute cluster, this simple
job scheduler can be used instead of LSF et al. It implements a single queue
of jobs, and it just tries to blindly execute them on the local CPU, without
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
    with 'MooseX::Daemonize';
    use Parallel::ForkManager;
    use Sys::CPU;
    use Cwd;
    use POSIX qw(ceil);
    
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::Schema;
    
    our $DEFAULT_CPUS = Sys::CPU::cpu_count();
    
    has 'deployment' => (is => 'ro',
                         isa => 'Str',
                         default => 'testing',
                         documentation => 'testing|production (default testing): determins which database is used to track locally scheduled jobs');
    
    has 'pidbase' => (is => 'rw',
                      isa => Dir,
                      coerce => 1,
                      builder => '_default_pidbase',
                      lazy => 1);
    
    has 'cpus' => (is => 'ro',
                   isa => 'Int',
                   default => $DEFAULT_CPUS);
    
    # for submit
    has 'o' => (is => 'ro',
                isa => 'Str',
                documentation => 'absolute path to stdout of scheduler (required for submit)');
    
    has 'e' => (is => 'ro',
                isa => 'Str',
                documentation => 'absolute path to stderr of scheduler (required for submit)');
    
    has 'a' => (is => 'ro',
                isa => 'Int',
                documentation => 'to specify a job array give the size of the array (default 1)',
                default => 1);
    
    
    sub BUILD {
        my $self = shift;
        my $d = $self->deployment;
        unless ($d eq 'testing' || $d eq 'production') {
            $self->throw("'$d' is not a valid deployment type; --deployment testing|production");
        }
        VRPipe::Persistent::SchemaBase->database_deployment($self->deployment);
    }
    
    method _default_pidbase {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler_output_root';
        return $vrp_config->$method_name();
    }
    
    after start {
        return unless $self->is_daemon;
        
        while (1) {
            $self->process_queue;
            sleep(5);
        }
    }
    
    method submit (Str $cmd) {
        my $o_file = $self->o;
        my $e_file = $self->e;
        unless ($o_file && $e_file) {
            $self->throw("-o and -e are required for submit");
        }
        
        my $array_size = $self->a;
        my $lsj = VRPipe::LocalSchedulerJob->create(cmd => $cmd, array_size => $array_size, cwd => cwd());
        
        my $user = getlogin || getpwuid($<);
        my @lsjs_args;
        foreach my $aid (1..$array_size) {
            my $this_o = $o_file;
            $this_o =~ s/\%I/$aid/g;
            my $this_e = $e_file;
            $this_e =~ s/\%I/$aid/g;
            push(@lsjs_args, { localschedulerjob => $lsj, aid => $aid, o_file => $this_o, e_file => $this_e, user => $user });
        }
        VRPipe::LocalSchedulerJobState->bulk_create_or_update(@lsjs_args);
        
        my $sid = $lsj->id;
        print "Job <$sid> is submitted\n";
    }
    
    method jobs (ArrayRef $ids) {
        my @lsjss;
        if (@$ids) {
            foreach my $id (@$ids) {
                push(@lsjss, $self->job_states($id));
            }
        }
        else {
            foreach my $lsj (VRPipe::LocalSchedulerJob->search({})) {
                push(@lsjss, $lsj->jobstates);
            }
        }
        
        return unless @lsjss;
        
        print join("\t", qw(JOBID INDEX STAT USER EXEC_HOST SUBMIT_TIME)), "\n";
        foreach my $lsjs (@lsjss) {
            my $lsj = $lsjs->localschedulerjob;
            print join("\t", $lsj->id, $lsjs->aid, $lsjs->current_status, $lsjs->user, $lsjs->host || '', $lsj->creation_time), "\n";
        }
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
            $self->throw("bad id format '$id'");
        }
        
        my ($lsj) = VRPipe::LocalSchedulerJob->search({id => $sid});
        
        unless ($lsj) {
            print "Job <$sid> is not found\n";
            return;
        }
        
        my @lsjss;
        foreach my $lsjs (VRPipe::LocalSchedulerJobState->search({localschedulerjob => $sid})) {
            if ($aid) {
                next unless $lsjs->aid == $aid;
            }
            push(@lsjss, $lsjs);
        }
        if ($aid && ! @lsjss) {
            print "Job <$sid\[$aid\]> is not found\n";
            return;
        }
        
        return @lsjss;
    }
    
    method unfinished_lsjss {
        return VRPipe::LocalSchedulerJobState->search({start_time => undef});
    }
    
    method process_queue {
        my @lsjss = $self->unfinished_lsjss;
        @lsjss || return;
        
        my $fm = Parallel::ForkManager->new($self->cpus);
        
        foreach my $lsjs (@lsjss) {
            $fm->start and next; # fork
            $lsjs->start_job;
            $fm->finish;
        }
        $fm->wait_all_children;
        
        return 1;
    }
    
    method kill (ArrayRef $ids) {
        my @lsjss;
        foreach my $id (@$ids) {
            push(@lsjss, $self->job_states($id));
        }
        
        foreach my $lsjss (@lsjss) {
            $lsjss->kill_job;
        }
    }
}