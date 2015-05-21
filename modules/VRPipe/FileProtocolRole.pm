
=head1 NAME

VRPipe::FileProtocolRole - a role that must be used by all FileProtocols

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

role VRPipe::FileProtocolRole {
    use Devel::GlobalDestruction;
    
    # the absolute path to the file (without any protocol prefixed) as a string
    has 'path' => (
        is       => 'ro',
        isa      => 'Str',
        required => 1
    );
    
    # the full protocol, including any username and password etc.
    has 'protocol' => (
        is       => 'ro',
        isa      => 'Str',
        required => 1
    );
    
    # returns our class name
    has 'pro' => (
        is      => 'ro',
        isa     => FileProtocol,
        lazy    => 1,
        builder => '_build_pro'
    );
    
    # implementations of open() should store the filehandle they generate in
    # here, allowing close() to work
    has _opened => (
        is  => 'rw',
        isa => 'Maybe[IO::File|FileHandle]'
    );
    
    method _build_pro {
        my $class = ref($self);
        my ($protocol) = $class =~ /.+::(.+)/;
        return $protocol;
    }
    
    # return a string that can be used in a command line to get the contents
    # of this file on to STDOUT for piping to something; it should do the right
    # thing if the file is compressed
    requires 'cat_cmd';
    
    # return a file handle on the opened file
    # first arg must be the open mode like '<' or '>'
    # second arg must be a hashref of possible options: permissions, backwards
    # options could be ignored if the protocol doesn't support them
    requires 'open';
    
    # below are some convenience methods that use required methods defined in
    # protocol subclasses
    
    method openr {
        return $self->open('<');
    }
    
    method openw {
        return $self->open('>');
    }
    
    method close {
        my $fh = $self->_opened || return;
        close($fh);
        $self->_opened(undef);
        return 1;
    }
    
    sub DEMOLISH {
        return if in_global_destruction;
        shift->close;
    }
}

1;
