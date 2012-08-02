
=head1 NAME

VRPipe::FileType::bam - bam filetype

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

bam files are the compressed binary forms of sam files:
L<http://samtools.sourceforge.net/> They hold aligned biological sequence data.

*** more documentation to come

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

class VRPipe::FileType::bam extends VRPipe::FileType::bin {
    our $correct_magic = [qw(037 213 010 004 000 000 000 000 000 377 006 000 102 103 002 000)];
    
    my $samtools_exe = file($ENV{SAMTOOLS}, 'samtools');
    
    around check_type {
        $self->$orig || return 0;
        #return $self->check_magic($self->file, $correct_magic); #*** for some reason, this, which in other contexts takes <1second, takes 10s of seconds in the Manager loop
        my $path = $self->file;
        return $path =~ /\.bam$/ ? 1 : 0;
    }
    
    method num_header_lines {
        my $path         = $self->file;
        my $headers      = `$samtools_exe view -H $path`;
        my @header_lines = split(/\n/, $headers);
        return scalar(@header_lines);
    }
    
    method num_records {
        my $path    = $self->file;
        my $records = `$samtools_exe view -c $path`;
        ($records) = $records =~ /^(\d+)/m;
        $records ||= 0;
        return $records;
    }
}

1;
