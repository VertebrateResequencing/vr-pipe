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
        isa     => 'Int',
        required => 1
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
        my $job_id = $self->job_id;
        my $job; # ... from db based on job_id ...
        my $cmd_line = 'echo wiggle'; # $job->cmd;
        return MooseX::Workers::Job->new(name => $job_id,
                                         command => sub { system("$cmd_line"); });
    }
    
    method _build_heartbeat_sub {
        return sub { print "heartbeat\n"; }
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
}

1;