=head1 NAME

VRPipe::DataSource::fofn - get pipeline input from a file of filenames

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

'fofn' stands for 'file of filenames', ie. the source should be a text file with
an absolute file path on each line. Each of these will become a
L<VRPipe::DataElement>.
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

class VRPipe::DataSource::fofn extends VRPipe::DataSource::list {
    method description {
        return "Use a simple list of absolute file paths (file-of-file-names) in a file as your source.";
    }
    method source_description {
        return "The path to a file with one absolute file path per line.";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will correspond to a single file from the file.";
        }
        
        return '';
    }
    
    method all (Defined :$handle!) {
        my @elements;
        foreach my $result ($self->_all_results(handle => $handle, skip_comments => 1, line_is_path => 1)) {
            my $withdraw = -e $result->{paths}->[0] ? 0 : 1;
            push(@elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => $result, withdrawn => $withdraw));
        }
        return \@elements;
    }
}

1;