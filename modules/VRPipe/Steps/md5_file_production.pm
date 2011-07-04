use VRPipe::Base;

class VRPipe::Steps::md5_file_production with VRPipe::StepRole {
    method inputs_definition {
        return { md5_file_input => VRPipe::StepIODefinition->get(type => 'any', max_files => -1, description => 'one or more files you want to calculate the md5 checksum of') };
    }
    method body_sub {
        return sub { my $self = shift;
                     my $input_files = $self->inputs->{md5_file_input};
                     @$input_files || return 1;
                     foreach my $vrfile (@$input_files) {
                        my $ifile = $vrfile->path;
                        my $ofile = $self->output_file(output_key => 'md5_files', basename => $ifile->basename.'.md5', type => 'txt', metadata => {md5_checksum_of => $ifile->stringify})->path;
                        $self->dispatch([qq{md5sum $ifile > $ofile}, $self->new_requirements(memory => 50, time => 1)]);
                     } };
    }
    method outputs_definition {
        return { md5_files => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1, description => '.md5 files for each input file containing its md5 checksum',
                                                            metadata => {md5_checksum_of => 'absoulte path of the file this .md5 file was generated for'}) };
    }
    method post_process_sub {
        return sub { my $self = shift;
                     my $output_files = $self->outputs->{md5_files};
                     foreach my $ofile (@$output_files) {
                        my $content = $ofile->slurp;
                        $content || return 0;
                        my ($md5) = split(" ", $content);
                        my $ifile = VRPipe::File->get(path => $ofile->metadata->{md5_checksum_of});
                        $ifile->md5($md5);
                        $ifile->update;
                     }
                     return 1; };
    }
    method description {
        return "Takes a file, calculates its md5 checksum, produces a file called <input filename>.md5, and updates the persistent database with the md5 of the file";
    }
}

1;
