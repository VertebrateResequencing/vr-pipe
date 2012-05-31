=head1 NAME

VRPipe::Persistent::Converter::mysql - a converter for MySQL

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

class VRPipe::Persistent::Converter::mysql with VRPipe::Persistent::ConverterRole { 

    method description {
        return "MySQL implemetation of database-specific methods used by Persistent";
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
			$datatype = 'text';
			$size = undef;
			return($datatype,$size,$is_numeric);
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
				$size = undef;
			}
		}
		return($datatype,$size,$is_numeric);
	}

    method get_datetime_type () {
		return 'datetime';
	}

    method get_boolean_type () {
		return 'bool';
	}

    method get_index_statements (Str $table_name, HashRef $for_indexing, Str $mode) {

		my @idx_cmds;
        my ($cols,$txt_cols);

		foreach my $k (keys %{$for_indexing}) {
			if ($for_indexing->{$k} eq 'text') {	# text index requires a size
				$txt_cols .= "$k(255),";
			}
			else {
				$cols .= "$k,";
			}
		}
		if ($cols) {
           if ($mode eq 'create') {
		       chop($cols);
		       push(@idx_cmds,"create index psuedo_idx on $table_name ($cols)");
           }
           else {
		       push(@idx_cmds,"drop index psuedo_idx on $table_name");
           }
		}
		if ($txt_cols) {
           if ($mode eq 'create') {
		       chop($txt_cols);
		       push(@idx_cmds,"create index txt_idx on $table_name ($txt_cols)");
           }
           else {
		       push(@idx_cmds,"drop index txt_idx on $table_name");
           }
		}
        return \@idx_cmds;
	}

	method get_index_cols (VRPipe::Persistent::Schema $schema, Str $table_name) {

		# Put columns in psuedo_idx or txt_idx into a hash

		use VRPipe::Config;
		my $vrp_config = VRPipe::Config->new();
		my $db_name = $vrp_config->production_dbname;

		my %idx_cols;
		my $rc  = $schema->storage->dbh_do(
			sub {
				my ($storage, $dbh, $table_name) = @_;
				my $res = $dbh->selectall_arrayref("show index from $table_name");
				foreach( @$res ) {
					if ($_->[2] eq 'txt_idx' or $_->[2] eq 'psuedo_idx') {
						my $col_name = $_->[4];
						my $res2 = $dbh->selectrow_hashref("select data_type from INFORMATION_SCHEMA.COLUMNS where TABLE_SCHEMA = '$db_name' and TABLE_NAME = '$table_name' and COLUMN_NAME = '$col_name'");
						$idx_cols{$col_name} = $res2->{data_type};
					}
				}
			},
			$table_name
		);
		return \%idx_cols;
	}

	method retype_index_cols (HashRef $idx_cols) {

		# Bools are implemented as tinyint; retype when checking if a table has changed during db upgrade

		foreach my $k (keys %{$idx_cols}) {
			if ($idx_cols->{$k} eq 'bool') {
				$idx_cols->{$k} = 'tinyint';
			}
		}
	}
}

1;
