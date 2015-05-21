
=head1 NAME

VRPipe::FileProtocol::file - handler for local files

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Handles local files with knowledge of how to get their file content streamed
out for piping purposes.

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

class VRPipe::FileProtocol::file with (VRPipe::FileProtocolRole, VRPipe::Base::FileMethods) {
    use Path::Class;
    use File::ReadBackwards;
    use IO::Uncompress::AnyUncompress;
    our $bgzip_magic = [37, 213, 10, 4, 0, 0, 0, 0, 0, 377, 6, 0, 102, 103, 2, 0];
    
    method cat_cmd {
        my $path = $self->path;
        if ($path =~ /\.gz$/) {
            return "zcat $path";
        }
        return "cat $path";
    }
    
    method open (OpenMode $mode, HashRef $args?) {
        my $path = file($self->path);
        my ($permissions, $backwards, $record_separator);
        if ($args) {
            $permissions      = $args->{permissions};
            $backwards        = $args->{backwards};
            $record_separator = $args->{record_separator};
        }
        
        $self->throw("Only modes <, > and >> are supported") unless $mode =~ /^(?:<|>)+$/;
        
        my $fh;
        
        # set up the open command, handling compressed files automatically
        my $open_cmd = $path;
        my $magic    = `file -bi $path`;
        ($magic) = split(';', $magic);
        if ($magic eq 'application/octet-stream' || $path =~ /\.gz$/) {
            if ($mode eq '<') {
                if ($backwards) {
                    $self->throw("Unable to read '$path' backwards when it is compressed");
                }
                
                # if it was made with Heng Li's bgzip it will be detected as a
                # gzip file, but will fail to be decompressed properly with
                # IO::Uncompress; manually detect the magic ourselves
                if ($self->check_magic($path, $bgzip_magic)) {
                    $open_cmd = "gunzip -c $path |";
                }
                else {
                    $fh = IO::Uncompress::AnyUncompress->new($path->stringify, AutoClose => 1);
                }
            }
            else {
                $open_cmd = "| gzip -c $mode $path";
            }
            
            open($fh, $open_cmd) unless $fh;
        }
        else {
            if ($mode eq '<' && $backwards) {
                # we'll open it with File::ReadBackwards, but first check it can
                # opened normally
                my $ok = open(my $testfh, '<', $path);
                if ($ok) {
                    close($testfh);
                    my @frb_args = ($path);
                    push(@frb_args, $record_separator) if $record_separator;
                    tie(*BW, 'File::ReadBackwards', @frb_args);
                    $fh = \*BW;
                }
            }
            else {
                my @args = ($mode);
                push(@args, $permissions) if $permissions;
                $fh = $path->open(@args);
            }
        }
        
        unless ($fh) {
            if ($mode eq '<') {
                $self->throw("Failed to open '$path': $!");
            }
            elsif ($mode eq '>') {
                $self->throw("Could not write to '$path': $!\n");
            }
        }
        
        $self->_opened($fh);
        return $fh;
    }
}

1;
