use VRPipe::Base;

class VRPipe::Steps::fastq_map with VRPipe::StepRole {
    method inputs_definition {
        return { md5_file_input => VRPipe::FileDefinition->get(name => 'file_input', type => 'any') };
    }
    method body_sub {
        return sub { my $self = shift;
                     my $input_files = $self->_input_files_from_input_or_element('md5_file_input');
                     @$input_files || return 1;
                     my $output_root = $self->output_root;
                     foreach my $ifile (@$input_files) {
                        my $ofile = Path::Class::File->new($output_root, $ifile->basename.'.md5');
                        $self->dispatch([qq{md5sum $ifile > $ofile}, $self->new_requirements(memory => 50, time => 1)]);
                     }
                     return 0; };
    }
    method post_process_sub {
        return sub { my $self = shift;
                     my $md5_files = $self->_input_files_from_input_or_element('md5_file_input');
                     my $output_root = $self->output_root;
                     foreach my $file (@$md5_files) {
                        my $ofile = Path::Class::File->new($output_root, $file->basename.'.md5');
                        my $content = $ofile->slurp;
                        $content || return 0;
                        my ($md5) = split(" ", $content);
                        my $vrfile = VRPipe::File->get(path => $file);
                        $vrfile->md5($md5);
                        $vrfile->update;
                     }
                     return 1; };
    }
    method outputs_definition {
        return { md5_files => VRPipe::FileDefinition->get(name => 'md5_files', type => 'txt', output_sub => sub { my ($self, $step) = @_;
                                                                                                                  my $input_files = $step->_input_files_from_input_or_element('md5_file_input');
                                                                                                                  my @md5_files;
                                                                                                                  foreach my $file (@$input_files) {
                                                                                                                      push(@md5_files, $file->basename.'.md5');
                                                                                                                  }
                                                                                                                  return [@md5_files]; }) };
    }
    method description {
        return "Takes a file, calculates its md5 checksum, produces a file called <input filename>.md5, and updates the persistent database with the md5 of the file";
    }
}

1;
