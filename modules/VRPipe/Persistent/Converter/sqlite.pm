=head1 NAME

VRPipe::Persistent::Converter::mysql - a converter for MySQL

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

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
	return 'datetime'; # datetime is magical; 'TEXT' and 'NUMERIC' do not work here, despite what the datatype docs imply
    }
    
    method get_boolean_type () {
	return 'INTEGER';
    }
    
    method get_index_statements (Str $table_name, HashRef $for_indexing, Str $mode) {
	my @idx_cmds;
	my $cols;
	
	foreach my $k (keys %{$for_indexing}) {
	    $cols .= "$k,";
	}
	if ($cols) {
	    if ($mode eq 'create') {
		chop($cols);
		push(@idx_cmds,"create index ${table_name}_psuedo_idx on $table_name ($cols)");
	    }
	    else {
		push(@idx_cmds,"drop index ${table_name}_psuedo_idx");
	    }
	}
	
	return \@idx_cmds;
    }
    
    method get_index_cols (VRPipe::Persistent::Schema $schema, Str $table_name) {
	# Put columns in psuedo_idx into a hash
	use VRPipe::Config;
	my $vrp_config = VRPipe::Config->new();
	my $db_name = $vrp_config->production_dbname;
	
	my %idx_cols;
	my $rc  = $schema->storage->dbh_do(
	    sub {
		my ($storage, $dbh, $table_name) = @_;
		my $res = $dbh->selectrow_hashref(qq[select sql from sqlite_master where tbl_name = '$table_name' and type = 'index' and name = '${table_name}_psuedo_idx']);
		my $create_index_sql = $res->{sql} || return; # not all classes have is_keys
		my ($cols) = $create_index_sql =~ /CREATE INDEX ${table_name}_psuedo_idx on $table_name \((.+)\)/i;
		foreach my $col_name (split(',', $cols)) {
		    my $res2 = $dbh->selectrow_hashref(qq[select sql from sqlite_master where name = '$table_name' and type = 'table']);
		    my $create_table_sql = $res2->{sql};
		    $create_table_sql =~ /$col_name (\w+)/;
		    $idx_cols{$col_name} = $1;
		}
	    },
	    $table_name
	);
	
	return \%idx_cols;
    }
}

1;
