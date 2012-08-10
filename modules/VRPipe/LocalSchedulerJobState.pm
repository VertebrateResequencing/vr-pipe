
=head1 NAME

VRPipe::LocalSchedulerJobState - job tracking for the local scheduler

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

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

class VRPipe::LocalSchedulerJobState extends VRPipe::Persistent {
    use Sys::Hostname;
    
    has 'localschedulerjob' => (is         => 'rw',
                                isa        => Persistent,
                                coerce     => 1,
                                belongs_to => 'VRPipe::LocalSchedulerJob',
                                traits     => ['VRPipe::Persistent::Attributes'],
                                is_key     => 1);
    
    has 'aid' => (is     => 'rw',
                  isa    => IntSQL [8],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_key => 1);
    
    has 'o_file' => (is     => 'rw',
                     isa    => File,
                     coerce => 1,
                     traits => ['VRPipe::Persistent::Attributes']);
    
    has 'e_file' => (is     => 'rw',
                     isa    => File,
                     coerce => 1,
                     traits => ['VRPipe::Persistent::Attributes']);
    
    has 'start_time' => (is          => 'rw',
                         isa         => Datetime,
                         coerce      => 1,
                         traits      => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    has 'end_time' => (is          => 'rw',
                       isa         => Datetime,
                       coerce      => 1,
                       traits      => ['VRPipe::Persistent::Attributes'],
                       is_nullable => 1);
    
    has 'exit_code' => (is          => 'rw',
                        isa         => IntSQL [5],
                        traits      => ['VRPipe::Persistent::Attributes'],
                        is_nullable => 1);
    
    has 'host' => (is          => 'rw',
                   isa         => Varchar [64],
                   traits      => ['VRPipe::Persistent::Attributes'],
                   is_nullable => 1);
    
    has 'user' => (is     => 'rw',
                   isa    => Varchar [64],
                   traits => ['VRPipe::Persistent::Attributes']);
    
    has 'pid' => (is          => 'rw',
                  isa         => IntSQL [6],
                  traits      => ['VRPipe::Persistent::Attributes'],
                  is_nullable => 1);
    
    __PACKAGE__->make_persistent();
    
    method current_status {
        if ($self->end_time) {
            if ($self->exit_code) {
                return 'EXIT';
            }
            else {
                return 'DONE';
            }
        }
        elsif ($self->start_time) {
            return 'RUN';
        }
        else {
            return 'PEND';
        }
    }
    
    method start_job {
        return if $self->start_time;
        
        $self->start_time(DateTime->now());
        $self->exit_code(undef);
        $self->update;
        
        my $cwd            = $self->localschedulerjob->cwd;
        my $submitters_env = $self->localschedulerjob->env;
        my $aid            = $self->aid;
        my $stdout_file    = $self->o_file;
        my $stderr_file    = $self->e_file;
        $stdout_file =~ s/\%I/$aid/g;
        $stderr_file =~ s/\%I/$aid/g;
        
        $self->disconnect;
        my $cmd_pid = fork();
        if (!defined $cmd_pid) {
            $self->throw("attempt to fork cmd failed: $!");
        }
        elsif ($cmd_pid) {
            # parent
        }
        elsif ($cmd_pid == 0) {
            # child
            $self->pid($$);
            $self->host(hostname());
            $self->user(getlogin || getpwuid($<));
            $self->update;
            
            # get all info from db and disconnect before using the info below
            my $cmd = $self->localschedulerjob->cmd;
            $stdout_file->dir->mkpath;
            my $sid        = $self->localschedulerjob->id;
            my $start_time = $self->start_time;
            my $host       = $self->host;
            my $user       = $self->user;
            $self->disconnect;
            
            chdir($cwd);
            open STDOUT, '>>', $stdout_file or $self->throw("Can't redirect STDOUT to '$stdout_file': $!");
            open STDERR, '>>', $stderr_file or $self->throw("Can't redirect STDERR to '$stderr_file': $!");
            
            print STDOUT "#--- vrpipe-local_scheduler report start ---\n# Job: ${sid}[$aid]\n# Started at: $start_time\n# Running on: $host for user: $user in working dir: $cwd\n# Cmd: $cmd\n# STDERR from this Cmd, if any, appears in $stderr_file\n# STDOUT from this Cmd, if any, appears below\n#---\n";
            
            while (my ($key, $val) = each %$submitters_env) {
                $ENV{$key} = $val;
            }
            $ENV{VRPIPE_LOCAL_JOBINDEX} = $aid;
            
            exec($cmd);
        }
        
        # wait for the cmd to finish
        my ($exit_code, $exit_meaning);
        if (waitpid($cmd_pid, 0) > 0) {
            $exit_code = $?;
            my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
            if ($core) {
                $exit_meaning = 'dumped core';
            }
            elsif ($sig == 9) {
                $exit_meaning = 'killed';
            }
            elsif ($rc) {
                my $message = $sig ? " after receiving signal $sig" : '';
                $exit_meaning = "exited with code $rc$message";
            }
            else {
                $exit_meaning = 'finished OK';
            }
        }
        else {
            $exit_code    = -1;
            $exit_meaning = 'lost track of PID';
        }
        
        # update db details for the job
        $self->end_time(DateTime->now());
        $self->exit_code($exit_code);
        $self->update;
        
        # write out info to stdout file
        open(my $ofh, '>>', $stdout_file) || $self->throw("Could not append to $stdout_file");
        my $end_time = $self->end_time;
        print $ofh "#---\n# Finished at: $end_time\n# Exit code: $exit_code ($exit_meaning)\n#--- vrpipe-local_scheduler report end ---\n\n";
        close($ofh);
    }
    
    method kill_job {
        return if $self->current_status eq 'DONE';
        
        my $pid = $self->pid;
        if ($pid) {
            if (!$self->end_time) {
                kill(9, $pid);
            }
        }
        
        $self->end_time(DateTime->now());
        $self->exit_code(22);
        $self->update;
    }
}
