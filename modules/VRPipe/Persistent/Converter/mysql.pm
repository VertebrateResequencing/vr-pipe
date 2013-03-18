
=head1 NAME

VRPipe::Persistent::Converter::mysql - a converter for MySQL

=head1 SYNOPSIS
        
        my $converter = VRPipe::Persistent::ConverterFactory->create('mysql', {});
        ($cname, $size, $is_numeric) = $converter->get_column_info(size => -1, is_numeric => 0);

=head1 DESCRIPTION

A Converter implementation for the MySQL database.

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

class VRPipe::Persistent::Converter::mysql with VRPipe::Persistent::ConverterRole {
    method description {
        return "MySQL implementation of database-specific methods used by Persistent";
    }
    
    method capitalise_name {
        return "MySQL";
    }
    
    method get_column_info (Int :$size!, Int :$is_numeric!) {
        my $datatype;
        
        # mysql has a limit on varchars of 255 up to version 5.0.3
        # size -1 indicates unlimited size
        
        if ($size > 255 || $size < 0) {
            $is_numeric = 0;
            $datatype   = 'text';
            $size       = undef;
            return ($datatype, $size, $is_numeric);
        }
        
        if ($is_numeric == 0) {
            $datatype = 'varchar';
        }
        else {
            if ($size < 3) {
                $datatype = 'tinyint';
            }
            elsif ($size < 5) {
                $datatype = 'smallint';
            }
            elsif ($size < 7) {
                $datatype = 'mediumint';
            }
            elsif ($size < 10) {
                $datatype = 'int';
            }
            else {
                $datatype = 'bigint';
                $size     = undef;
            }
        }
        
        return ($datatype, $size, $is_numeric);
    }
    
    method get_datetime_type () {
        return 'datetime';
    }
    
    method get_boolean_type () {
        return 'tinyint';
    }
    
    method get_index_statements (Str $table_name, HashRef $for_indexing, Str $mode) {
        my @idx_cmds;
        foreach my $col (sort keys %{$for_indexing}) {
            my $index_name = $table_name . '_idx_' . $col;
            
            my $spec = $col;
            if ($for_indexing->{$col} eq 'text') {
                # text index requires a size
                $spec = "$col(255)";
            }
            
            if ($mode eq 'create') {
                push(@idx_cmds, "create index $index_name on $table_name ($spec)");
            }
            else {
                push(@idx_cmds, "drop index $index_name on $table_name");
            }
        }
        return \@idx_cmds;
    }
    
    method get_index_cols (VRPipe::Persistent::Schema $schema, Str $table_name) {
        # Put columns in psuedo_idx or txt_idx into a hash
        use VRPipe::Config;
        my $vrp_config = VRPipe::Config->new();
        my $db_name    = $vrp_config->production_dbname;
        
        my %idx_cols;
        my $rc = $schema->storage->dbh_do(
            sub {
                my ($storage, $dbh, $table_name) = @_;
                
                # first check we even have this table
                my $res = $dbh->selectall_arrayref("show tables like '$table_name'");
                return if @$res == 0;
                
                $res = $dbh->selectall_arrayref("show index from $table_name");
                foreach (@$res) {
                    if ($_->[2] =~ /${table_name}_idx_/) {
                        my $col_name = $_->[4];
                        my $res2     = $dbh->selectrow_hashref("select data_type from INFORMATION_SCHEMA.COLUMNS where TABLE_SCHEMA = '$db_name' and TABLE_NAME = '$table_name' and COLUMN_NAME = '$col_name'");
                        $idx_cols{$col_name} = $res2->{data_type};
                    }
                }
            },
            $table_name
        );
        
        return \%idx_cols;
    }
    
    method index_creation_style {
        return 'single';
    }
    
    method get_isolation_read_committed_sql {
        return 'SET TRANSACTION ISOLATION LEVEL READ COMMITTED';
    }
}

1;
