
=head1 NAME

VRPipe::FileList - store an unordered list of VRPipe::File objects

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

You probably want to always use get() to safely create or get FileLists.

Direct use of create() will always create and return a new FileList, even if
one with the same set of files already exists.

These lists are immuatble.

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

class VRPipe::FileList extends VRPipe::Persistent {
    __PACKAGE__->make_persistent(has_many => [members => 'VRPipe::FileListMember']);
    
    our $empty_fl;
    
    around get (ClassName|Object $self: Persistent :$id?, ArrayRef[VRPipe::File] :$files?) {
        $self->throw("You cannot supply both id and files") if $id && $files;
        
        if ($id) {
            return $self->$orig(id => $id);
        }
        elsif ($files) {
            if (@$files == 0) {
                # we need to return the same keyvallist every time, but the
                # below code doesn't find those with no members, so we need to
                # special-case
                unless ($empty_fl) {
                    #*** I have no idea how to do this as single sql query in DBIx::Class...
                    my $pager = VRPipe::FileList->get_column_values_paged('id', {});
                    PAGES: while (my $ids = $pager->next) {
                        foreach my $fl_id (@$ids) {
                            my $members = VRPipe::FileListMember->search({ filelist => $fl_id }, { rows => 1 });
                            unless ($members) {
                                $empty_fl = VRPipe::FileList->get(id => $fl_id);
                                last PAGES;
                            }
                        }
                    }
                }
                return $empty_fl if $empty_fl;
            }
            
            my %flids;
            my $index = 0;
            foreach my $file (@$files) {
                ++$index;
                foreach my $flid (VRPipe::FileListMember->get_column_values('filelist', { file => $file->id })) {
                    $flids{$flid}++;
                }
            }
            
            my $expected_count = @$files;
            foreach my $flid (sort { $a <=> $b } keys %flids) {
                next unless $flids{$flid} == $expected_count;
                next unless VRPipe::FileListMember->search({ filelist => $flid }) == $expected_count;
                return $self->$orig(id => $flid);
            }
            return $self->create(files => $files);
        }
        else {
            $self->throw("FileList->get requires id or files");
        }
    }
    
    around create (ClassName|Object $self: ArrayRef[VRPipe::File] :$files!) {
        # create a new row, then use the new id to create new
        # FileListMember rows for each supplied file
        my $list = $self->$orig();
        
        my $index = 0;
        my @flm_args;
        foreach my $file (@$files) {
            push(@flm_args, { filelist => $list->id, file => $file->id });
        }
        VRPipe::FileListMember->bulk_create_or_update(@flm_args);
        
        return $list;
    }
    
    method files {
        my @return;
        foreach my $member ($self->members) {
            push(@return, $member->file);
        }
        return @return;
    }
}

1;
