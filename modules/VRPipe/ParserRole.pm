
=head1 NAME

VRPipe::ParserRole - a role required for parsers

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Parsers are based around the idea of getting a single record (eg. a single
line, depending on the file format) from a file and parsing it in fields and
making those fields available in a parsed record data structure (usually an
array reference).

The reference to the data structure does not change (only its contents), so you
can efficiently call next_record and just keep accessing the same data
structure in a loop.

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

role VRPipe::ParserRole {
    has 'file' => (
        is       => 'ro',
        isa      => VRPFileOrHandle,
        coerce   => 1,
        required => 1
    );
    
    has 'type' => (
        is      => 'ro',
        isa     => ParserType,
        lazy    => 1,
        builder => '_build_type'
    );
    
    has 'parsed_record' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        default => sub { [] }
    );
    
    has 'fh' => (
        is      => 'ro',
        isa     => AnyFileHandle,
        writer  => '_set_fh',
        lazy    => 1,
        builder => '_get_fh'
    );
    
    has 'filename' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_get_filename'
    );
    
    has '_buffer_store' => (
        is      => 'ro',
        traits  => ['Array'],
        isa     => 'ArrayRef[Str]',
        default => sub { [] },
        handles => {
            '_pushback'     => 'push',
            '_getback'      => 'shift',
            '_empty_buffer' => 'clear'
        }
    );
    
    has '_vrpipe_file' => (
        is     => 'ro',
        isa    => 'VRPipe::File',
        writer => '_set_vrpipe_file'
    );
    
    has '_tell' => (
        is  => 'rw',
        isa => 'Int'
    );
    
    has '_current_results' => (
        is  => 'rw',
        isa => 'ArrayRef|HashRef'
    );
    
    has '_header_parsed' => (
        is        => 'ro',
        isa       => 'Int',
        writer    => '_set_header_parsed',
        predicate => '_got_header'
    );
    
    method BUILD {
        $self->_get_header;
    }
    
    method _build_type {
        my $class = ref($self);
        my ($type) = $class =~ /.+::(.+)/;
        return $type;
    }
    
    method _get_fh {
        my $file = $self->file;
        if (ref($file) && (File->check($file) || $file->isa('VRPipe::File'))) {
            my $vrpf = $file->isa('VRPipe::File') ? $file : VRPipe::File->create(path => $file->absolute, FileType->check($self->type) ? (type => $self->type) : ());
            $self->_set_vrpipe_file($vrpf); # because on destruction, $vrpf will close the opened filehandle
            $vrpf->disconnect;
            return $vrpf->openr();
        }
        else {
            return $file;
        }
    }
    
    method _get_filename {
        $self->fh;
        my $vrpf = $self->_vrpipe_file;
        if ($vrpf) {
            my $filename = $vrpf->path->stringify;
            $vrpf->disconnect;
            return $filename;
        }
        else {
            return '-';
        }
    }
    
    requires 'next_record';
    
    method _readline (AnyFileHandle $fh) {
        my $line = $self->_getback;
        unless (defined $line) {
            $line = <$fh>;
        }
        return $line;
    }

=head2 _save_position
 
 Title   : _save_position
 Usage   : $self->_save_position()
 Function: Internal method for parser authors. Saves the current filehandle
           position; for use before seeking.
 Returns : boolean (true if this is a seekable filehandle)
 Args    : n/a

=cut
    
    method _save_position {
        my $fh = $self->fh() || return;
        
        my $tell = tell($fh);
        if ($tell == -1) {
            $self->warn("this parsing method doesn't work on piped input");
            return;
        }
        my $current_results = ref($self->parsed_record) eq 'ARRAY' ? [@{ $self->parsed_record }] : { %{ $self->parsed_record } };
        
        $self->_tell($tell);
        $self->_current_results($current_results);
        
        return 1;
    }

=head2 _get_header
 
 Title   : _get_header
 Usage   : $self->_get_header()
 Function: Internal method for parser authors. This does nothing. If your file
           format has a header, you should implement it by overriding this.
           In your implementation, _set_header_parsed() should be called after
           successfully parsing your header.
           It should return true if the header was parsed, false if not.
 Returns : boolean
 Args    : n/a

=cut
    
    method _get_header {
        unless ($self->_header_parsed) {
            $self->_set_header_parsed(0);
        }
        return 0;
    }

=head2 _set_header_parsed
 
 Title   : _set_header_parsed
 Usage   : $self->_set_header_parsed()
 Function: Internal method for parser authors. Ask if the header of the current
           file has been parsed.
 Returns : boolean
 Args    : n/a (optionally, supply a filehandle position as from tell() that
           denotes the start of the first record after the header; by default
           it assumes the current filehandle position matches this criteria)

=cut
    
    around _set_header_parsed (Int $tell?) {
        my $fh = $self->fh() || return;
        unless (defined $tell) {
            $tell = tell($fh);
        }
        $self->$orig($tell);
    }

=head2 seek
 
 Title   : seek
 Usage   : $self->seek($pos, $whence);
 Function: Behaves exactly like Perl's standard seek(), except that if the
           filehandle was made by opening a .gz file for reading, you can
           effectively seek backwards.
 Returns : boolean (for success)
 Args    : position to seek to, position to seek from

=cut
    
    method seek (Int $tell, Int $whence) {
        my $fh = $self->fh() || return;
        
        if (ref($fh) eq 'IO::Uncompress::Gunzip') {
            # we can't go backwards, so close and re-open without changing our
            # fh id
            close($fh);
            my $z = IO::Uncompress::Gunzip->new($self->file->stringify);
            $self->{_fh} = $z;
            $fh = $z;
            $self->_set_fh($fh);
        }
        
        CORE::seek($fh, $tell, $whence);
    }

=head2 _seek_first_record
 
 Title   : _seek_first_record
 Usage   : $self->_seek_first_record()
 Function: Internal method for parser authors. Seeks back to before the first
           result (ie. after a header if present) so that next_result() will
           behave as if it was called for the first time on a file.
           You should call _save_position() first, then call this, then do your
           work, then call _restore_position().
 Returns : n/a
 Args    : n/a

=cut
    
    method _seek_first_record {
        my $tell;
        if ($self->_got_header) {
            $tell = $self->_header_parsed;
            if ($tell == -1) {
                $self->warn("this parsing method doesn't work on piped input");
                return;
            }
        }
        else {
            $tell = 0;
        }
        
        $self->seek($tell, 0);
    }

=head2 _restore_position
 
 Title   : _restore_position
 Usage   : $self->_restore_position()
 Function: Internal method for parser authors. Restores the current filehandle
           position to the location when _save_position() was last called; for
           use after you've seeked somewhere and done some work.
           If _save_position() was called at the true start of the file, this
           actually calls _seek_first_record().
 Returns : n/a
 Args    : n/a

=cut
    
    method _restore_position {
        my $fh = $self->fh() || return;
        
        if ($self->_tell == 0) {
            # we might have saved position before parsing the header, but now have
            # the flag set that we've parsed the header; don't seek back before
            # the header!
            $self->_seek_first_record();
            return;
        }
        
        $self->seek($self->_tell, 0);
        if (ref($self->parsed_record) eq 'ARRAY') {
            my @current_results = @{ $self->_current_results };
            for my $i (0 .. $#current_results) {
                $self->parsed_record->[$i] = $current_results[$i];
            }
        }
        else {
            while (my ($key, $val) = each %{ $self->_current_results }) {
                $self->parsed_record->{$key} = $val;
            }
        }
        
        return 1;
    }
}

1;
