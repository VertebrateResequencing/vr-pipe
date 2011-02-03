=head1 NAME

VRPipe::Base::Cleanable - instance cleanup methods

=head1 SYNOPSIS

package Foo;
use Moose;


=head1 DESCRIPTION

Methods to do things when an instance of a VRPipe object is destroyed.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

role VRPipe::Base::Cleanable {
    has 'cleanup_methods'  => ( is => 'ro',
                                isa => 'HashRef[Str]',
                                default => sub { {} },
                                auto_deref => 1,
                                writer => '_set_cleanup_methods' );
    has 'files_to_unlink'  => ( is => 'ro',
                                isa => 'HashRef[Str]',
                                default => sub { {} },
                                auto_deref => 1,
                                writer => '_set_files_to_unlink' );
    
=head2 register_for_cleanup

 Title   : register_for_cleanup
 Usage   : $self->register_for_cleanup($method_name);
 Function: Store a method that will be called upon object destruction.
 Returns : n/a
 Args    : method name string

=cut
    method register_for_cleanup (Str $method!) {
        if ($self->can($method)) {
            my %cleanup_methods = $self->cleanup_methods;
            $cleanup_methods{$method} = 1;
            $self->_set_cleanup_methods(\%cleanup_methods);
        }
    }
    
=head2 unregister_for_cleanup

 Title   : unregister_for_cleanup
 Usage   : $self->unregister_for_cleanup($method_name);
 Function: Unstore a method previously registered so that it won't be called
           upon object destruction.
 Returns : n/a
 Args    : method name string

=cut
    method unregister_for_cleanup (Str $method!) {
        my %cleanup_methods = $self->cleanup_methods;
        delete $cleanup_methods{$method};
        $self->_set_cleanup_methods(\%cleanup_methods);
    }
    
=head2 register_for_unlinking

 Title   : register_for_unlinking
 Usage   : $self->register_for_unlinking($filename);
 Function: Store a filepath that will be unlinked upon object destruction.
 Returns : n/a
 Args    : filepath(s)

=cut
    method register_for_unlinking (@files) {
        foreach my $file (@files) {
            my %files_to_unlink = $self->files_to_unlink;
            $files_to_unlink{$file} = 1;
            $self->_set_files_to_unlink(\%files_to_unlink);
        }
    }
    
=head2 unregister_for_unlinking

 Title   : unregister_for_unlinking
 Usage   : $self->unregister_for_unlinking($filename);
 Function: Unstore a file previously registered so that it won't be unlinked
           upon object destruction.
 Returns : n/a
 Args    : filepath(s)

=cut
    method unregister_for_unlinking (@files) {
        foreach my $file (@files) {
            my %files_to_unlink = $self->files_to_unlink;
            delete $files_to_unlink{$file};
            $self->_set_files_to_unlink(\%files_to_unlink);
        }
    }
    
=head2 DEMOLISH

 Title   : DEMOLISH
 Usage   : n/a
 Function: Our DEMOLISH method calls all our cleanup methods and unlinks
           registered files.
 Returns : n/a
 Args    : n/a

=cut
    method DEMOLISH {
        my %cleanup_methods = $self->cleanup_methods;
        foreach my $method (keys %cleanup_methods) {
            $self->$method();
        }
        
        my %files_to_unlink = $self->files_to_unlink;
        foreach my $file (keys %files_to_unlink) {
            unlink($file);
        }
    }
    
=head2 Signal handling

    By default, DESTROY is called for each object when __DIE__, SIGINT, or
    SIGTERM is received. In case that the default DESTROY handlers are not
    flexible enough, do NOT override $SIG{INT} or $SIG{TERM} handlers. Either
    create a child of Base and modify the DESTROY handler, or add static
    set_die_handler and unset_die_handler calls to Base. 
    
    The signals should be probably trapped only when requested, not by default.
    Maybe the subroutine register_for_unlinking should do this, as it this was
    introduced only to really unlink certain files?

=cut
    our $SIGNAL_CAUGHT_EVENT = "Signal caught.\n";
    $SIG{TERM} = sub { die $SIGNAL_CAUGHT_EVENT; };
    $SIG{INT}  = sub { die $SIGNAL_CAUGHT_EVENT; };
}

1;
