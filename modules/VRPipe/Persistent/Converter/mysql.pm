use VRPipe::Base;

class VRPipe::Persistent::Converter::mysql with VRPipe::Persistent::ConverterRole
{ 

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

    method get_index_statements (Str $table_name, HashRef $for_indexing ) {

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
		   chop($cols);
		   push(@idx_cmds,"create index psuedo_idx on $table_name ($cols)");
		}
		if ($txt_cols) {
		   chop($txt_cols);
		   push(@idx_cmds,"create index txt_idx on $table_name ($txt_cols)");
		}
        return \@idx_cmds;
	}

}

1;
