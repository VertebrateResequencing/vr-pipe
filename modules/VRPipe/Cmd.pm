=head1 NAME

VRPipe::Cmd - execute a job in a subprocess with a heartbeat process in the bg

=head1 SYNOPSIS

use VRPipe::Cmd;

my $cmd = VRPipe::Cmd->new(job_id => 2445);
$cmd->run;

# or

perl -M VRPipe::Cmd -e 'VRPipe::Cmd->new->run' --job_id 2445


=head1 DESCRIPTION

This takes the id of a Job (supplied to new() or via command line args) and
executes its associated command line, outputting stdout and stderr to the
special VRPipe locations, and running a heart-beat process in the background so
that you can check in the database it hasn't silently stalled.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

class VRPipe::Cmd with VRPipe::Base::LivingProcesses {
    use MooseX::Workers::Job;
    use Time::Format;
    
    has job_id => (
        is      => 'ro',
        isa     => Persistent,
        required => 1
    );
    
    has job => (
        is      => 'ro',
        isa     => 'VRPipe::Job',
        lazy => 1,
        builder => '_build_job'
    );
    
    has stdout => (
        is      => 'ro',
        isa     => FileOrHandle,
        builder => '_build_stdout',
        coerce  => 1,
        lazy    => 1,
        required => 1
    );
    
    has _stdofh => (
        is      => 'rw',
        isa     => 'Maybe[FileHandle]'
    );
    
    has stderr => (
        is      => 'ro',
        isa     => FileOrHandle,
        builder => '_build_stderr',
        coerce  => 1,
        lazy    => 1,
        required => 1
    );
    
    has _stdefh => (
        is      => 'rw',
        isa     => 'Maybe[FileHandle]'
    );
    
    method _build_stdout {
        return '~/.vrpipe.stdout';
    }
    method _build_stderr {
        return '~/.vrpipe.stderr';
    }
    
    method _std_fh (Str $type) {
        my $fh_method_name = "_std${type}fh";
        my $fh = $self->$fh_method_name();
        
        unless ($fh) {
            my $file_method_name = $type eq 'o' ? 'stdout' : 'stderr';
            my $file = $self->$file_method_name();
            if (ref($file) eq 'Path::Class::File') {
                $fh = $file->open('>>');
                if (-s $file) {
                    print $fh "\n\n";
                }
                print $fh "---\{$time{'yyyy/mm/dd'} $time{'hh:mm:ss'}\}---\n\n";
            }
            else {
                $fh = $file;
            }
            $self->$fh_method_name($fh);
        }
        
        return $fh;
    }
    method _stdout_fh {
        return $self->_std_fh('o');
    }
    method _stderr_fh {
        return $self->_std_fh('e');
    }
    
    method _build_processes {
        my $job = $self->job;
        my $cmd_line = $job->cmd;
        warn "setting main mwjob to a $job with cmd [$cmd_line]\n";
        return MooseX::Workers::Job->new(name => $self->job_id,
                                         command => sub { my $exit_code = system($cmd_line); exit($exit_code); });
    }
    
    method _build_job {
        my $job_id = $self->job_id;
        eval "require VRPipe::Job;";
        return VRPipe::Job->get(id => $job_id);
    }
    
    method _build_heartbeat_sub {
        return sub { my $self = shift || return; my $job = $self->job; $job->heartbeat(DateTime->now()); $job->update; }
    }
    
    method worker_stdout (Str $output, MooseX::Workers::Job $job) {
        my $fh = $self->_stdout_fh;
        chomp($output);
        print $fh $output, "\n";
    }
    method worker_stderr (Str $output, MooseX::Workers::Job $job) {
        my $fh = $self->_stderr_fh;
        chomp($output);
        print $fh $output, "\n";
    }
    
    method worker_started (MooseX::Workers::Job $mwjob) {
        warn sprintf("%s(%s,%s) started\n", $mwjob->name, $mwjob->ID, $mwjob->PID);
        return if $self->_is_heartbeat_worker($mwjob);
        warn "this was the main worker\n";
    }
    
    method worker_done (MooseX::Workers::Job $mwjob) {
        warn sprintf("%s(%s,%s) finished\n", $mwjob->name, $mwjob->ID, $mwjob->PID);
        return if $self->_is_heartbeat_worker($mwjob);
        warn "this was the main worker\n";
    }
    
    method worker_error (MooseX::Workers::Job $mwjob) {
        warn sprintf("%s(%s,%s) had an error\n", $mwjob->name, $mwjob->ID, $mwjob->PID);
        return if $self->_is_heartbeat_worker($mwjob);
        warn "this was the main worker\n";
    }
    
    method worker_manager_stop {
        if ($self->_stdofh) {
            $self->_stdofh->close;
            $self->_stdofh(undef);
        }
        if ($self->_stdefh) {
            $self->_stdefh->close;
            $self->_stdefh(undef);
        }
    }
    
    method sig_TERM {
        warn "got sig_TERM(@_)\n";
        if ($self->num_workers > 0) {
            my $engine = $self->Engine;
            my @ids = sort { $a <=> $b } $engine->get_worker_ids;
            foreach my $id (@ids) {
                my $worker = $engine->get_worker($id) || next;
                warn "will try and kill worker $id\n";
                $worker->kill();
            }
        }
    }
    
    method DEMOLISH {
        warn "Cmd being demolished\n";
        if ($self->num_workers > 0) {
            my $engine = $self->Engine;
            my @ids = sort { $a <=> $b } $engine->get_worker_ids;
            foreach my $id (@ids) {
                my $worker = $engine->get_worker($id) || next;
                warn "will try and kill worker $id\n";
                $worker->kill();
            }
        }
    }
}

1;