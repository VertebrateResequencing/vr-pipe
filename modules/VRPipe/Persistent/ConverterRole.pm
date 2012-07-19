
=head1 NAME

VRPipe::Persistent::ConverterRole - a role required of all SQL converters

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

role VRPipe::Persistent::ConverterRole {
    requires 'description';
    requires 'capitalise_name';
    requires 'get_column_info';
    requires 'get_datetime_type';
    requires 'get_boolean_type';
    requires 'get_index_statements';
    requires 'get_index_cols';
    requires 'index_creation_style';
    
    method index_statements (VRPipe::Persistent::Schema $schema, Str $mode) {
        my @idx_cmds;
        foreach my $class (keys %{ $schema->class_mappings }) {
            my $table_name = $class;
            $table_name =~ s/.*:://;
            $table_name = lc($table_name);
            
            my $meta         = $class->meta;
            my $for_indexing = $meta->get_attribute('cols_to_idx')->get_value($meta);
            
            if (keys %{$for_indexing}) {
                push(@idx_cmds, @{ $self->get_index_statements($table_name, $for_indexing, $mode) });
            }
        }
        return \@idx_cmds;
    }
}

1;
