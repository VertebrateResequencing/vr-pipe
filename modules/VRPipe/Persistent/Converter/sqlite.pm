
=head1 NAME

VRPipe::Persistent::Converter::mysql - a converter for SQLite

=head1 SYNOPSIS
        
        my $converter = VRPipe::Persistent::ConverterFactory->create('sqlite', {});
        ($cname, $size, $is_numeric) = $converter->get_column_info(size => -1, is_numeric => 0);

=head1 DESCRIPTION

A Converter implementation for the SQLite database.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Persistent::Converter::sqlite with VRPipe::Persistent::ConverterRole {
    method description {
        return "SQlite implemetation of database-specific methods used by Persistent";
    }
    
    method capitalise_name {
        return "sqlite";
    }
    
    method get_column_info (Int :$size!, Int :$is_numeric!) {
        my $datatype;
        if ($is_numeric == 0) {
            $datatype = 'TEXT';
        }
        else {
            $datatype = 'INTEGER';
        }
        return ($datatype, undef, $is_numeric);
    }
    
    method get_datetime_type () {
        return 'datetime';  # datetime is magical; 'TEXT' and 'NUMERIC' do not work here, despite what the datatype docs imply
    }
    
    method get_boolean_type () {
        return 'INTEGER';
    }
    
    method get_index_statements (Str $table_name, HashRef $for_indexing, Str $mode) {
        my @idx_cmds;
        foreach my $col (sort keys %{$for_indexing}) {
            my $index_name = $table_name . '_idx_' . $col;
            
            if ($mode eq 'create') {
                push(@idx_cmds, "create index $index_name on $table_name ($col)");
            }
            else {
                push(@idx_cmds, "drop index $index_name");
            }
        }
        return \@idx_cmds;
    }
    
    method get_index_cols (VRPipe::Persistent::Schema $schema, Str $table_name) {
        # Put columns in psuedo_idx into a hash
        use VRPipe::Config;
        my $vrp_config = VRPipe::Config->new();
        my $db_name    = $vrp_config->production_dbname;
        
        my %idx_cols;
        my $rc = $schema->storage->dbh_do(
            sub {
                my ($storage, $dbh, $table_name) = @_;
                my $res = $dbh->selectall_arrayref(qq[select sql from sqlite_master where tbl_name = '$table_name' and type = 'index' and name like '${table_name}_idx_%']);
                foreach (@$res) {
                    my ($col_name) = $_->[0] =~ /CREATE INDEX ${table_name}_idx_\w+ on $table_name \((.+)\)/i;
                    
                    my $res2             = $dbh->selectrow_hashref(qq[select sql from sqlite_master where name = '$table_name' and type = 'table']);
                    my $create_table_sql = $res2->{sql};
                    $create_table_sql =~ /$col_name (\w+)/;
                    $idx_cols{$col_name} = $1;
                }
            },
            $table_name
        );
        
        return \%idx_cols;
    }
    
    method index_creation_style {
        return 'all';
    }
    
    method get_isolation_change_sql (Bool :$repeatable_read = 0) {
        return if $repeatable_read; # by default it is serializable, which is an acceptable substitute
        # else, we want read committed...
        #*** hmmm, can SQLite do read-committed-type behaviour? If not, can we
        # support it?
        return;
    }
}

1;
