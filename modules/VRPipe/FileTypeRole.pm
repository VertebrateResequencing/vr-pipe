
=head1 NAME

VRPipe::FileTypeRole - a role that must be used by all FileTypes

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

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

role VRPipe::FileTypeRole {
    use VRPipe::Parser;
    
    has 'file' => (
        is       => 'rw',
        isa      => File,
        coerce   => 1,
        required => 1
    );
    
    has 'type' => (
        is      => 'ro',
        isa     => FileType,
        lazy    => 1,
        builder => '_build_type'
    );
    
    has 'record_separator' => (
        is      => 'ro',
        isa     => 'Maybe[Str]',
        builder => '_build_record_separator'
    );
    
    has 'read_backwards' => (
        is      => 'ro',
        isa     => 'Bool',
        builder => '_build_read_backwards'
    );
    
    has 'parser' => (
        is      => 'ro',
        does    => 'VRPipe::ParserRole',
        lazy    => 1,
        builder => '_build_parser'
    );
    
    method _build_type {
        my $class = ref($self);
        my ($type) = $class =~ /.+::(.+)/;
        return $type;
    }
    
    method _build_record_separator {
        return;
    }
    
    method _build_read_backwards {
        return 0;
    }
    
    method _build_parser {
        return VRPipe::Parser->create($self->type, { file => $self->file });
    }
    
    method check_type {
        my $file = $self->file;
        my $type = $self->type;
        $file =~ s/\.gz$// unless $type eq 'gz';
        if ($file =~ /\.$type$/) {
            return 1;
        }
        return 0;
    }
    
    method num_lines {
        return $self->num_header_lines + $self->num_records;
    }
    
    method num_header_lines {
        return 0;
    }
    
    method num_records {
        my $path = $self->file;
        my $cat = $path =~ /\.gz$/ ? 'zcat' : 'cat';
        open(my $wc, "$cat $path | wc -l |") || $self->throw("$cat $path | wc -l did not work");
        my ($lines) = split(" ", <$wc>);
        CORE::close($wc);
        return $lines;
    }
}

1;
