
=head1 NAME

VRPipe::Base::Debuggable - debug message handling

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Provides warning, error, debug and log message handling. *** more documentation
to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011,2013 Genome Research Limited.

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

role VRPipe::Base::Debuggable {
    use Carp qw(carp cluck confess);
    use Time::Format;
    use Cwd qw(getcwd);
    use File::Spec;
    use IO::Capture::Stderr;
    
    # MooseX::ClassAttribute won't let us set things globally across all classes
    # that compose this role, and class_has attributes don't get picked up by
    # MooseX::StrictConstructor, so we don't use it.
    our $GLOBAL_VERBOSITY;
    our $GLOBAL_WRITE_LOGS = 0;
    our $GLOBAL_LOG_FILE = File::Spec->catfile($ENV{HOME}, '.VRPipe.log');
    
    # these are just to allow verbose and log_file to be set during new()
    has '_verbose' => (
        is       => 'ro',
        isa      => VerbosityValue,
        default  => 0,
        init_arg => 'verbose',
        trigger  => \&verbose
    );
    has '_log_file' => (
        is       => 'ro',
        isa      => 'Str|Path::Class::File',
        init_arg => 'log_file',
        trigger  => \&log_file
    );

=head2 verbose
 
 Title   : verbose
 Usage   : $obj->verbose(1); # per-instance
           VRPipe::Any::Class->set_verbose_global(2); # global override -
                                           # make all instances of all classes
                                           # behave like they have the given
                                           # verbosity
           VRPipe::Any::Class->clear_verbose_global; # disables global override
 Function: Sets verbose level for how ->warn() and ->debug() behave
           -1 = no warning, no debug messages
            0 = standard small warning, no debug messages
            0.5 = tiny warning with no line numbers
            1 = warning with stack trace
            2 = warning becomes throw
 Returns : The current verbosity setting (number between -1 and 2)
 Args    : none to get, a verbosity value to set: -1,0,0.5,1 or 2

=cut
    
    method verbose (ClassName|Object $self: VerbosityValue $value?, VerbosityValue $old_val?) {
        unless (ref($self)) {
            return $GLOBAL_VERBOSITY || 0;
        }
        
        if (defined $value) {
            $self->{verbose} = ($value);
        }
        
        if (defined $GLOBAL_VERBOSITY) {
            return $GLOBAL_VERBOSITY;
        }
        
        return $self->{verbose} || 0;
    }
    
    method set_verbose_global (ClassName|Object $self: VerbosityValue $value) {
        $GLOBAL_VERBOSITY = $value;
    }
    
    method clear_verbose_global (ClassName $class:) {
        undef $GLOBAL_VERBOSITY;
    }

=head2 warn
 
 Title   : warn
 Usage   : $obj->warn("warning message");
 Function: Warns about something; strength of warning determined by ->verbose().
           If logging is turned on will also output warning to log file.
 Returns : n/a
 Args    : A string giving a warning message

=cut
    
    method warn (ClassName|Object $self: Str $message) {
        my $verbose = $self->verbose();
        return if $verbose <= -1;
        
        if ($verbose >= 2) {
            $self->throw($message);
        }
        elsif ($verbose == 1) {
            cluck($message);
        }
        elsif ($verbose == 0.5) {
            CORE::warn($message . "\n");
        }
        else {
            carp($message);
        }
        
        $self->log($message);
    }

=head2 debug
 
 Title   : debug
 Usage   : $obj->debug("This is debugging output");
 Function: Warns a debugging message when verbose is > 0.
           If logging is turned on will also output message to log file.
 Returns : n/a
 Args    : Message string to debug about

=cut
    
    method debug (ClassName|Object $self: Str $message) {
        if ($self->verbose > 0) {
            $self->log($message);
            $message .= "\n" unless $message =~ /\n$/;
            CORE::warn "$time{'hh:mm:ss'} $$ | ", $message;
        }
    }

=head2 throw
 
 Title   : throw
 Usage   : $obj->throw("throwing exception message");
 Function: Throws an exception, which, if not caught with an eval or
           a try block will provide a nice stack trace to STDERR
           with the message. (Uses Carp's confess().) Automatically includes
           the date and current working directory.
           If logging is turned on will also output throw message to log file.
 Returns : n/a
 Args    : A string giving a descriptive error message

=cut
    
    method throw (ClassName|Object $self: Str $message = '[no message]') {
        my $cwd = getcwd();
        #my $first_line = "FATAL ERROR on $time{'yyyy/mm/dd hh:mm:ss'} in $cwd";
        my $t          = time();
        my $first_line = "FATAL ERROR at epoch $t in $cwd";
        
        my $throw_message = "\n-----------------------------------------------------------------------------\n";
        $throw_message .= $first_line . "\n\n" . $message . "\n\n" . $self->stack_trace;
        $throw_message .= "\n-----------------------------------------------------------------------------\n\n";
        
        $self->log($throw_message, 2);
        
        die $throw_message;
    }
    
    sub stack_trace {
        my $capture = IO::Capture::Stderr->new();
        $capture->start();
        cluck();
        $capture->stop();
        
        my @confess;
        foreach my $read ($capture->read) {
            foreach my $line (split(/\n/, $read)) {
                $line =~ s/\t/    /g;
                
                if ($line =~ /^    \S+Debuggable::throw.+ called at (.+)/) {
                    $line = "Thrown from $1";
                }
                elsif ($line =~ /^ at \S+Debuggable.pm line \d+/ || $line =~ /^    \S+Debuggable::stack_trace.+ called at/) {
                    next;
                }
                
                push(@confess, $line);
            }
        }
        
        return join("\n", @confess);
    }

=head2 write_logs
 
 Title   : write_logs
 Usage   : $obj->write_logs(1);
 Function: Turn the writing of log messages on or off. Changing this acts
           globally across all VRPipe instances of any class.
 Returns : boolean
 Args    : none to get, boolean to set

=cut
    
    method write_logs (ClassName|Object $self: Bool $set?) {
        if (defined $set) {
            $GLOBAL_WRITE_LOGS = $set;
        }
        return $GLOBAL_WRITE_LOGS;
    }

=head2 log_file
 
 Title   : log_file
 Usage   : $obj->log_file("/path/to/logfile");
 Function: Get or set the path to the log file. Changing this acts globally
           across all VRPipe instances of any class.
 Returns : string
 Args    : none to get, string to set

=cut
    
    method log_file (ClassName|Object $self: Str|Path::Class::File $set?, Str|File $old_val?) {
        if (defined $set) {
            open(my $fh, ">>", $set) || $self->throw("Could not append to log file '$set'");
            close($fh);
            $GLOBAL_LOG_FILE = $set;
        }
        return $GLOBAL_LOG_FILE;
    }

=head2 log
 
 Title   : log
 Usage   : $obj->log('message');
 Function: Appends message to file set in log_file() if write_logs() is on.
 Returns : n/a
 Args    : log message, optional int to decide the strength of the message, as
           per verbose(), using current value of verbose() as default.

=cut
    
    method log (ClassName|Object $self: Str $message, VerbosityValue $verbose?) {
        return unless $self->write_logs();
        $verbose ||= $self->verbose;
        return unless $verbose > -1;
        
        my $prefix = '';
        
        $self->{_log_count}++;
        if ($self->{_log_count} == 1) {
            $prefix = "----------\n$time{'yyyy/mm/dd'}\n----------\n";
        }
        
        $prefix .= "$time{'hh:mm:ss'} $$ | ";
        
        # for some reason Carp's shortmess and longmess methods aren't returning
        # the same thing as what carp()/croak() produce! hacky work-around...
        # (we die instead of throw since throw calls log()...)
        my $capture = IO::Capture::Stderr->new();
        $capture->start();
        $verbose >= 1 ? cluck($message) : carp($message);
        $capture->stop();
        
        my $file = $self->log_file;
        my $opened = open(my $fh, ">>", $file);
        if ($opened) {
            print $fh $prefix, $capture->read, "\n";
            close($fh);
        }
        else {
            carp("Unable to write to log file '$file', disabling logging!");
            $self->write_logs(0);
        }
    }
}

1;
