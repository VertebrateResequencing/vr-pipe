=head1 NAME

VRPipe::Base::FileMethods - commonly used methods for dealing with files

=head1 SYNOPSIS

use MooseX::Declare;
class VRPipe::MyClass with (VRPipe::Base::FileMethods) {
    #...
    
    method my_method (Str $input, Str $source) {
        my $path = $self->catfile($self->cwd, $input);
        my $ofh = $self->open('>', $path);
        my $ifh = $self->open('<', "/path/to/input.txt.gz");
        $self->copy($source, $dest);
    }
}

=head1 DESCRIPTION

Provides a grab-bag of very commonly used file-related methods from external
(CPAN or CORE) modules, so you can call $self->common_method instead of { use
external::module; external::module->common_method }. Also adds some wrapping
to make error handling consistent for all the methods.

Allows for the possibilty to reimplement these methods here, to reduce the
number of external dependencies. Though currently all methods come from CORE
modules (in recent versions of Perl at least), so no need.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

role VRPipe::Base::FileMethods with (VRPipe::Base::Debuggable) {
    use MooseX::Aliases;
    use VRPipe::Base::Types qw(FileNameOrHandle);
    use File::Spec;
    use Cwd qw(abs_path cwd);
    use File::Basename;
    use File::Path qw(make_path remove_tree);
    use File::Copy;
    
=head2 cwd

 Title   : cwd
 Usage   : my $path = $obj->cwd(); 
 Function: Get the current working directory. Just an alias to Cwd::cwd.
 Returns : string
 Args    : n/a

=cut
    method cwd () {
        return Cwd::cwd;
    }
    
=head2 make_path

 Title   : make_path (alias mkpath)
 Usage   : $obj->make_path($path);
 Function: Make directories, like mkdir -p. An alias to File::Path::make_path,
           but with automatic VRPipe-style handling of verbosity and errors.
 Returns : n/a
 Args    : as per File::Path::make_path

=cut
    method make_path (Str $path, @args) {
        my $args;
        if (@args && ref($args[-1])) {
            $args = pop @args;
        }
        else {
            $args = {};
        }
        
        unless (defined $args->{verbose}) {
            $args->{verbose} = $self->verbose;
        }
        $args->{error} = \my $err;
        
        push(@args, $args);
        
        File::Path::make_path($path, @args);
        
        if (@$err) {
            my $messages = '';
            for my $diag (@$err) {
                my ($file, $message) = %$diag;
                if ($file eq '') {
                    $messages .= "make_path general error: $message\n";
                }
                else {
                    $messages .= "make_path problem with $file: $message\n";
                }
            }
            $self->throw($messages);
        }
    }
    alias mkpath => 'make_path';
    
=head2 remove_tree

 Title   : remove_tree (alias rmtree)
 Usage   : $obj->remove_tree($path);
 Function: Remove a directory structe, like rm -rf. An alias to
           File::Path::remove_tree, but with automatic VRPipe-style handling of
           verbosity and errors.
 Returns : n/a
 Args    : as per File::Path::make_path

=cut
    method remove_tree (Str $path, @args) {
        my $args;
        if (@args && ref($args[-1])) {
            $args = pop @args;
        }
        else {
            $args = {};
        }
        
        unless (defined $args->{verbose}) {
            $args->{verbose} = $self->verbose;
        }
        $args->{error} = \my $err;
        
        push(@args, $args);
        
        File::Path::remove_tree($path, @args);
        
        if (@$err) {
            my $messages = '';
            for my $diag (@$err) {
                my ($file, $message) = %$diag;
                if ($file eq '') {
                    $messages .= "remove_tree general error: $message\n";
                }
                else {
                    $messages .= "remove_tree problem with $file: $message\n";
                }
            }
            $self->throw($messages);
        }
    }
    alias rmtree => 'remove_tree';
    
=head2 copy

 Title   : copy (alias cp)
 Usage   : $obj->copy($source, $dest);
 Function: Copy a file. An alias to File::Copy::copy, but with VRPipe-style
           handling of errors. Does not return a boolean; throws on failure
           instead.
 Returns : n/a
 Args    : as per File::Copy::copy

=cut
    method copy (FileOrHandle $source, FileOrHandle $dest) {
        my $success = File::Copy::copy($source, $dest);
        unless ($success) {
            $self->throw("copy of $source => $dest failed: $!");
        }
    }
    alias cp => 'copy';
    
=head2 move

 Title   : move (alias mv)
 Usage   : $obj->move($source, $dest);
 Function: Move a file. An alias to File::Copy::move, but with VRPipe-style
           handling of errors. Does not return a boolean; throws on failure
           instead.
 Returns : n/a
 Args    : as per File::Copy::move

=cut
    method move (FileOrHandle $source, FileOrHandle $dest) {
        my $success = File::Copy::move($source, $dest);
        unless ($success) {
            $self->throw("move of $source => $dest failed: $!");
        }
    }
    alias mv => 'move';
}

1;
