
=head1 NAME

VRPipe::DataSource::list - get pipeline input from a text file

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This DataSource will be of limited use. Each resulting L<VRPipe::DataElement>
will consist of a line from the source text file, but they will not be
considered paths; most pipelines need file paths to work with. *** more
documentation to come

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

class VRPipe::DataSource::list with VRPipe::DataSourceTextRole {
    method description {
        return "Use a simple list of items in a file as your source.";
    }
    
    method source_description {
        return "The path to a file with one item per line.";
    }
    
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will correspond to a single line from the file.";
        }
        
        return '';
    }
    
    method _open_source {
        my $file = $self->source_file;
        return $file->openr;
    }
    
    method all (Defined :$handle!, Bool :$skip_comments = 1, Bool :$line_is_path = 0) {
        my @element_args;
        foreach my $result ($self->_all_results(handle => $handle, skip_comments => $skip_comments, line_is_path => $line_is_path)) {
            push(@element_args, { datasource => $self->_datasource_id, result => $result });
        }
        $self->_create_elements(\@element_args);
    }
    
    method _all_results (Defined :$handle!, Bool :$skip_comments = 1, Bool :$line_is_path = 0) {
        my $key_name = $line_is_path ? 'paths' : 'line';
        
        my @results;
        while (<$handle>) {
            if (/^\s*$/) {
                next;
            }
            elsif ($skip_comments && /^#/) {
                next;
            }
            chomp;
            
            my $result = $_;
            $result = $line_is_path ? [VRPipe::File->create(path => file($result)->absolute)->path->stringify] : $result; # we can't bulk_create VRPipe::Files because they do fancy stuff duing create()
            push(@results, { $key_name => $result });
        }
        
        return @results;
    }
}

1;
