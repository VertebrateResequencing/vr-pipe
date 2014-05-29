
=head1 NAME

VRPipe::Base::FileMethods - commonly used methods for dealing with files

=head1 SYNOPSIS
    
    use VRPipe::Base;
    
    class VRPipe::MyClass with (VRPipe::Base::FileMethods) {
        #...
        
        method my_method (Str $input, Str $source, Str $dest) {
            my $tempdir = $self->tempdir();
            ($handle, $tempfile) = $self->tempfile();
        }
    }

=head1 DESCRIPTION

Provides a grab-bag of very commonly used file-related methods from external
(CPAN or CORE) modules, so you can call $self->common_method instead of { use
external::module; external::module->common_method }. Also adds some wrapping to
make error handling consistent for all the methods.

Allows for the possibilty to reimplement these methods here, to reduce the
number of external dependencies.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

role VRPipe::Base::FileMethods {
    use MooseX::Aliases;
    use Cwd qw(cwd);
    use File::Path qw(make_path remove_tree);
    use File::Copy;
    use File::Temp;
    use Digest::MD5;
    
    our $cat_marker = "---------------------------------VRPipe--concat---------------------------------\n";
    
    # File::Temp will auto-delete temporary files and dirs when our instance is
    # destroyed, but we need to keep a reference to them until then to stop
    # them being deleted too early
    has _file_temps => (
        traits  => ['Array'],
        is      => 'ro',
        isa     => 'ArrayRef[File::Temp::Dir|File::Temp::File]',
        lazy    => 1,
        default => sub { [] },
        handles => { _remember_file_temp => 'push' }
    );

=head2 cwd
 
 Title   : cwd
 Usage   : my $path = $obj->cwd(); 
 Function: Get the current working directory.
 Returns : string
 Args    : n/a

=cut
    
    method cwd () {
        return Dir(Cwd::cwd);
    }

=head2 make_path
 
 Title   : make_path (alias mkpath)
 Usage   : $obj->make_path($path);
 Function: Make directories, like mkdir -p. An alias to File::Path::make_path,
           but with automatic VRPipe-style handling of verbosity and errors.
 Returns : n/a
 Args    : as per File::Path::make_path

=cut
    
    method make_path (Dir $path, @args) {
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
    
    method remove_tree (Dir $path, @args) {
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
        
        #*** we need to update_stats_from_disc for everything in the db under
        #    this path
        
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

=head2 tempfile
 
 Title   : tempfile
 Usage   : my ($handle, $tempfile) = $obj->tempfile(); 
 Function: Get a temporary filename and a handle opened for writing and
           and reading. Just an alias to File::Temp::tempfile.
 Returns : a list consisting of temporary handle and temporary filename
 Args    : as per File::Temp::tempfile

=cut
    
    method tempfile {
        my $ft = File::Temp->new(@_);
        $self->_remember_file_temp($ft);
        return ($ft, Path::Class::File->new($ft->filename));
    }
    alias tmpfile => 'tempfile';

=head2 tempdir
 
 Title   : tempdir
 Usage   : my $tempdir = $obj->tempdir(); 
 Function: Creates and returns the name of a new temporary directory. Just an
           alias to File::Temp::newdir that defaults to storing within the
           logging directory.
 Returns : The name of a new temporary directory.
 Args    : as per File::Temp::newdir

=cut
    
    method tempdir {
        shift;
        my ($template, %args);
        if (@_ >= 3) {
            ($template, %args) = @_;
        }
        elsif (@_ == 2) {
            %args = @_;
        }
        elsif (@_ == 1) {
            $template = shift;
        }
        unless (defined $args{DIR} || $args{TMPDIR}) {
            eval "require VRPipe::Persistent::SchemaBase;";
            eval "require VRPipe::Config;";
            my $vrp_config  = VRPipe::Config->new();
            my $method_name = VRPipe::Persistent::SchemaBase->database_deployment . '_logging_directory';
            $args{DIR} = Path::Class::Dir->new($vrp_config->$method_name(), 'tmp');
            unless (-d $args{DIR}) {
                $self->make_path($args{DIR});
            }
        }
        
        my $ft = File::Temp->newdir($template ? ($template) : (), %args);
        $self->_remember_file_temp($ft);
        return Path::Class::Dir->new($ft->dirname);
    }
    alias tmpdir => 'tempdir';

=head2 concatenate
 
 Title   : concatenate
 Usage   : $obj->concatenate($source, $destination,
                             unlink_source => 1,
                             add_marker => 1); 
 Function: append the content of $source file to the end of $destination file
           with a marker (as understood by the VRPipe::Parser::cat parser)
           appended at the end. Optionally delete the source file afterwards.
 Returns : n/a
 Args    : VRPipe::File $source, VRPipe::File $destination, optionally
           unlink_source => Bool (default false), add_marker => Bool (default
           true) and max_lines => Int (no max by default; if set, upto half this
           value lines will be copied over from the start of the source, and
           upto half this value from the end of the source)

=cut
    
    method concatenate (VRPipe::File $source, VRPipe::File $destination, Bool :$unlink_source = 0, PositiveInt :$max_lines?, Bool :$add_marker = 1) {
        # it's important that we have up-to-date e and s for source and dest
        $source->update_stats_from_disc;
        $destination->update_stats_from_disc;
        
        my $copy_over = $source->s;
        if ($destination->s) {
            if ($add_marker) {
                my $last_line = $destination->last_line;
                unless ($last_line && $last_line eq $cat_marker) {
                    $self->add_cat_marker($destination);
                }
            }
        }
        elsif ($unlink_source) {
            $source->move($destination) if $source->e;
            $copy_over     = 0;
            $unlink_source = 0;
        }
        
        if ($copy_over) {
            my $ifh               = $source->open('<', backwards => 0);
            my $ofh               = $destination->open('>>');
            my $copied_lines      = 0;
            my $first_half_limit  = int($max_lines / 2) if $max_lines;
            my $second_half_limit = $max_lines - $first_half_limit if $max_lines;
            
            my $do_tail = 0;
            while (<$ifh>) {
                if ($first_half_limit && ($copied_lines >= $first_half_limit)) {
                    $do_tail++;
                    last if $do_tail > $second_half_limit;
                }
                else {
                    print $ofh $_;
                    $copied_lines++;
                }
            }
            
            if ($do_tail) {
                print $ofh "--[lines skipped during VRPipe-concat]--\n" unless $do_tail == $second_half_limit;
                close($ifh);
                my $path = $source->path;
                open($ifh, "tail -n $second_half_limit $path |") || $self->throw("Could not tail $path\n");
                while (<$ifh>) {
                    print $ofh $_;
                }
                close($ifh);
            }
            
            # *** could do with checks on expected line numbers...
            
            $source->close;
            $destination->close;
        }
        
        $self->add_cat_marker($destination) if $add_marker;
        
        if ($unlink_source) {
            $source->remove;
        }
    }
    
    method add_cat_marker (VRPipe::File $file) {
        my $fh = $file->open('>>');
        print $fh $cat_marker;
        $file->close;
    }
    
    method check_magic (File $path, ArrayRef[Int] $correct_magic) {
        # *** hardly cross-platform?!
        
        # *** od can hang for thousands of seconds, so at least disconnect from
        #     the db before doing this
        #*** how?! we are a base class with no db handling... instead all our
        # callers would have to do the disconnect...
        # Hopefully using -N 16 will avoid the hangs in future?...
        
        my $magic = `od -N 16 -b $path`;
        
        my (undef, @magic) = split(/\s/, $magic);
        my $is_correct = 1;
        foreach my $m (@$correct_magic) {
            my $this_m = shift(@magic);
            if ($this_m != $m) {
                $is_correct = 0;
                last;
            }
        }
        return $is_correct;
    }
    
    method verify_md5 (File $path, Str $md5) {
        my $vrfile = VRPipe::File->create(path => $path);
        
        if ($self->file_md5($vrfile) eq $md5) {
            $vrfile->md5($md5);
            $vrfile->update;
            return 1;
        }
        else {
            return 0;
        }
    }
    
    method file_md5 (VRPipe::File $vrfile) {
        my $path = $vrfile->path;
        $self->throw("Cannot get the md5 of a file ($path) that does not exist") unless -e $path;
        $vrfile->disconnect;
        open(my $fh, '<', $path) || $self->throw("Could not open file $path: $!");
        binmode($fh);
        return Digest::MD5->new->addfile($fh)->hexdigest;
    }
    
    method hashed_dirs (Str $hashing_string, PositiveInt $levels = 4) {
        my $dmd5 = Digest::MD5->new();
        $dmd5->add($hashing_string);
        my $md5 = $dmd5->hexdigest;
        my @chars = split("", $md5);
        return @chars[0 .. $levels - 1];
    }
}

1;
