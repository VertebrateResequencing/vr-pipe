=head1 NAME

VRPipe::DataSourceTextRole - a role for DataSources that work with text files

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

role VRPipe::DataSourceTextRole with VRPipe::DataSourceRole {
    has 'source_file' => => (is => 'ro',
                             isa => 'VRPipe::File',
                             lazy => 1,
                             builder => '_build_source_file');
    
    method _build_source_file {
        my $source = file($self->source)->absolute;
        return VRPipe::File->create(path => $source, type => 'txt');
    }
    
    method _has_changed {
        my $old_md5 = $self->_changed_marker || return 1;
        my $file = $self->source_file;
        my $new_md5 = $self->file_md5($file);
        if ($new_md5 ne $old_md5) {
            $file->update_md5($new_md5);
            return 1;
        }
        else {
            return 0;
        }
    }
    
    method _update_changed_marker {
        my $file = $self->source_file;
        my $md5 = $self->file_md5($file);
        $file->update_md5($md5);
        $self->_changed_marker($md5);
    }
}

1;