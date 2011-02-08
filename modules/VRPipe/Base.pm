=head1 NAME

VRPipe::Base - a convenience role for all VRPipe modules

=head1 SYNOPSIS

use MooseX::Declare;

class VRPipe::MyClass with VRPipe::Base {
    #...
    
    method my_method (Str $input) {
        # throw an exception
        $self->throw("This is an exception") if $bad_thing_happened;

        # catch and handle a possible exception
        eval {
            my $other = VRPipe::othermodule->new()
            $other->dangerous_method();
        };
        if ($@) {
            # got an exception, do something about it...
        }
        else {
            # didn't get an exception, do something else...
        }

        # warn about something
        $self->warn("This is a warning") if $something_not_quite_right;

        # print a debugging message for those who want to read it
        $self->debug("This is a debugging message");

        # change the verbosity of another object
        my $other = VRPipe::othermodule->new(verbose => 1);
        $other->verbose(2);

        # turn on logging of the above kinds of messages
        $self->write_logs(1);
        $self->log_file($self->catfile(qw(my special logfile location)));

        # only log a message to the log file; don't print to STDERR
        $self->log('my log message');

        # if you do something that will require special consideration when
        # destroying ourselves, write a method to handle that and register it
        $self->register_for_cleanup('my_cleanup_method');
    }
}

=head1 DESCRIPTION

Provides common needs for many VRPipe modules: ties together debug-related
message handling, temp file cleanup and commonly needed CPAN methods.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Base with (VRPipe::Base::Debuggable, VRPipe::Base::Cleanable, VRPipe::Base::FileMethods) {
    use MooseX::StrictConstructor;
}

1;
