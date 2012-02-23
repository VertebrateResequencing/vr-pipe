use VRPipe::Base;

role VRPipe::Persistent::ConverterRole {

    requires 'description';
    requires 'capitalise_name';
    requires 'get_column_info';
    requires 'get_datetime_type';
    requires 'get_boolean_type';
    requires 'get_index_statements';

    method index_statements (VRPipe::Persistent::Schema $schema) {

		my @idx_cmds;

        foreach my $class (keys %{$schema->class_mappings}) {
            my $table_name = $class;
            $table_name =~ s/.*:://;
            $table_name = lc($table_name);

            my $meta = $class->meta;
            my $for_indexing = $meta->get_attribute('idx_keys')->get_value($meta);

            if (keys %{$for_indexing}) {
				push(@idx_cmds, @{$self->get_index_statements($table_name, $for_indexing)});
            }
        }
        return \@idx_cmds;
    }

}

1;
