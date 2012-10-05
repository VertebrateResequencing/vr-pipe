
=head1 NAME

VRPipe::Runner - execute a command and know if it is running or not

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A Runner is used to run some command. It lets you check if the command is
running, or if at least it ought to run soon (if it was scheduled to run by a
job management system). It's state tracking beyond that is deliberately simple.
It does not remember if the command you want run already completed successfully
before, nor does it track when and if your command completes, in failure or
otherwise.

You're supposed to use it for commands you simply want to keep running. These
commands should also not care about what directory they're run from (ie. any
paths should be absolute).

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Runner extends VRPipe::Persistent::Living {
    has 'cmd' => (
        is     => 'rw',
        isa    => Text,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    __PACKAGE__->make_persistent();
    
    method running {
        if ($self->alive && $self->heartbeat) {
            return 1;
        }
        return 0;
    }
    
    method pending {
        if ($self->alive && !$self->heartbeat) {
            return 1;
        }
        return 0;
    }
    
    method run {
        # start our heartbeat so that another process can check if we're
        # running, and also make sure we end if murdered by another process
        $self->start_beating;
        $self->react_to_being_murdered;
        
        # now that we've done our job, erase our existence
        $self->commit_suicide;
    }
    
    # {
    #     my $cmd_pid = fork();
    #     my $heartbeat_pid;
    #     if (!defined $cmd_pid) {
    #         $self->throw("attempt to fork cmd failed: $!");
    #     }
    #     elsif ($cmd_pid) {
    #         # parent, start a heartbeat process to monitor the cmd
    #         $heartbeat_pid = fork();
    #         if (!defined $heartbeat_pid) {
    #             $self->throw("attempt to fork for heartbeat failed: $!");
    #         }
    #         elsif ($heartbeat_pid == 0) {
    #             # child, initiate a heartbeat that will end when the cmd stops
    #             # running
    #             my $interval = $self->heartbeat_interval;
    #             sleep(1);
    #             while (1) {
    #                 my $still_running = kill(0, $cmd_pid);
    #                 $self->heartbeat(DateTime->now());
    #                 $self->update;
    #                 $self->disconnect;
    #                 sleep $interval;
    #             }
    #             exit(0);
    #         }
    #     }
    #     elsif ($cmd_pid == 0) {
    #         # child, actually run the command after changing to the correct dir
    #         # and redirecting output to files
    #         $self->pid($$);
    #         $self->host(hostname());
    #         $self->user(getlogin || getpwuid($<));
    #         $self->update;
    
    #         # get all info from db and disconnect before using the info below
    #         my $dir         = $self->dir;
    #         my $cmd         = $self->cmd;
    #         my $stdout_file = $self->stdout_file->path; #*** call here in child, then below in parent, creating 2 identical rows in db?... should not be possible!
    #         $stdout_file->dir->mkpath;
    #         my $stderr_file = $self->stderr_file->path; #*** but only 1 copy of the e?!
    #         $self->disconnect;
    
    #         open STDOUT, '>', $stdout_file or $self->throw("Can't redirect STDOUT to '$stdout_file': $!");
    #         open STDERR, '>', $stderr_file or $self->throw("Can't redirect STDERR to '$stderr_file': $!");
    #         chdir($dir);
    
    #         # exec is supposed to get our $cmd to run whilst keeping the same
    #         # $cmd_pid, but on some systems like Ubuntu the sh (dash) is a bit
    #         # sucky: http://www.perlmonks.org/?node_id=785284
    #         # We can't force list mode in the normal way because we actually
    #         # require the use of the shell to do things like run multi-line
    #         # commands and pipes etc. $cmd having a different pid to $cmd_pid
    #         # matters because we need to know the correct pid if the server
    #         # needs to kill it later. Instead we force the use of a better shell
    #         # (default bash) for everything, which might be less efficient in
    #         # some cases, but the difference is going to be meaningless for us.
    #         my $shell = VRPipe::Config->new->exec_shell;
    #         if ($shell) {
    #             exec {$shell} $shell, '-c', $cmd;
    #         }
    #         else {
    #             exec $cmd;
    #         }
    #     }
    
    #     # wait for the cmd to finish
    #     my $exit_code = $self->_wait_for_child($cmd_pid);
    #     $self->reselect_values_from_db;
    
    #     # update db details for the job
    #     $self->end_time(DateTime->now());
    #     $self->stdout_file->update_stats_from_disc(retries => 3);
    #     $self->stderr_file->update_stats_from_disc(retries => 3);
    #     my $o_files = $self->output_files;
    #     if (@$o_files) {
    #         foreach my $o_file (@$o_files) {
    #             $o_file->update_stats_from_disc(retries => 3);
    #         }
    #     }
    #     elsif ($stepstate) {
    #         $stepstate->update_output_file_stats;
    #     }
    #     $self->running(0);
    #     $self->finished(1);
    #     $self->exit_code($exit_code);
    #     $self->update;
    
    #     # reap the heartbeat
    #     kill(9, $heartbeat_pid);
    #     $self->_wait_for_child($heartbeat_pid);
    # }
    
    # method _wait_for_child (Int $pid) {
    #     if (waitpid($pid, 0) > 0) {
    #         return $?;
    #         #*** below code should be used elsewhere by some kind of front-end
    #         #    reporting method, which would access the raw exit code from
    #         #    ->exit_code instead of $?
    #         #my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
    #         #if ($core) {
    #         #    $self->warn("$identifier dumped core");
    #         #    return $core;
    #         #}
    #         #elsif ($sig == 9) {
    #         #    $self->warn("$identifier was killed");
    #         #    return $sig;
    #         #}
    #         #elsif ($rc) {
    #         #    my $message = $sig ? " after receiving signal $sig" : '';
    #         #    $self->warn("$identifier exited with code $rc$message");
    #         #    return $rc;
    #         #}
    #         #else {
    #         #    #warn "child $pid finished OK ($?)\n";
    #         #    return $rc;
    #         #}
    #     }
    #     else {
    #         my $identifier = "child process $pid for command [" . $self->cmd . "]";
    #         $self->warn("$identifier disappeared!");
    #     }
    # }
}

1;
