
=head1 NAME

VRPipe::FileList - store an unordered list of VRPipe::File objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

You probably want to always use get() to safely create or get FileLists.

Direct use of create() will always create and return a new FileList, even if
one with the same set of files already exists.

These lists are immutable.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::FileList extends VRPipe::Persistent with VRPipe::PersistentListRole {
    sub _member_class { 'VRPipe::FileListMember' }
    sub _member_key   { 'file' }
    sub _foreign_key  { 'filelist' }
    
    __PACKAGE__->make_persistent(has_many => [members => 'VRPipe::FileListMember']);
    
    around get (ClassName|Object $self: Persistent :$id?, ArrayRef[VRPipe::File] :$files?) {
        return $self->_get_list($orig, $id, $files);
    }
    
    around create (ClassName|Object $self: ArrayRef[VRPipe::File] :$files!) {
        return $self->_create_list($orig, $files);
    }
    
    method files {
        return $self->_instantiated_members;
    }
}

1;
