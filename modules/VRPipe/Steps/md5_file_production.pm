use VRPipe::Base;

class VRPipe::Steps::md5_file_production with VRPipe::StepRole {
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
                     my $md5_files = $self->outputs->{md5_files};
                     foreach my $vrfile (@$md5_files) {
                        my $content = $vrfile->slurp;
                        $content || return 0;
                        my ($md5) = split(" ", $content);
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

=pod

my $pipeline1_output_def = VRPipe::FileDefinition->get(name => 'element_output_file',
                                                       type => 'txt',
                                                       output_sub => sub  { my ($self, undef, undef, $data_element) = @_;
                                                                            return $data_element->result.'.o'; });

my $single_step = VRPipe::Step->get(name => 'element_outputter',
                                    inputs_definition => { },
                                    body_sub => sub { my $self = shift;
                                                      my $element_name = $self->data_element->result;
                                                      my $ofile = $self->outputs->{the_only_output}->path;
                                                      $self->dispatch([qq{sleep 5; echo "output for $element_name" > $ofile}, $self->new_requirements(memory => 60, time => 1)]);
                                                      return 0; },
                                    post_process_sub => sub { return 1 },
                                    outputs_definition => { the_only_output => $pipeline1_output_def },
                                    description => 'outputs the data element result to a file');
                                    
=cut