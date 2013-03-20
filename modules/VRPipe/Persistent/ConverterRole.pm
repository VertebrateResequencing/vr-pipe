
=head1 NAME

VRPipe::Persistent::ConverterRole - a role required of all SQL converters

=head1 SYNOPSIS
        
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {}); # eg mysql
        
        $idx_sql = $converter->index_statements($self, $mode);
        
        ($cname, $size, $is_numeric) = $converter->get_column_info(size => -1, is_numeric => 0);
        $cname = $converter->get_boolean_type();

=head1 DESCRIPTION

Converters are required for each supported database, and implement various
methods to allow Persistent objects to be created, and the Schema to be
deployed, in a database independent way.

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
    requires 'get_isolation_change_sql';
    
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
    
    method get_db_schema_version (VRPipe::Persistent::Schema $schema) {
        my $db_version;
        my $rc = $schema->storage->dbh_do(
            sub {
                my ($storage, $dbh) = @_;
                my $res = $dbh->selectall_arrayref("select version from dbix_class_deploymenthandler_versions order by id desc");
                $db_version = $res->[0][0];
            }
        );
        
        return $db_version;
    }
}

1;
