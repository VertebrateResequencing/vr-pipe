use VRPipe::Base;

class VRPipe::Steps::fastq_split with VRPipe::StepRole {
    method inputs_definition {
        return { fastq_files => VRPipe::FileDefinition->get(name => 'fastq_files', type => 'fastq') };
    }
    method body_sub {
        return sub { my $self = shift;
                     my $input_files = $self->inputs->{fastq_files};
                     @$input_files || return 1;
                     my $output_root = $self->output_root;
                     foreach my $vrfile (@$input_files) {
                        my $ifile = $vrfile->path;
                        my $ofile = Path::Class::File->new($output_root, $ifile->basename.'.md5');
                        VRPipe::File->get(path => $ofile, type => "fastq");
                        $self->dispatch([qq{md5sum $ifile > $ofile}, $self->new_requirements(memory => 50, time => 1)]);
                     } };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method outputs_definition {
        return { fastq_files => VRPipe::FileDefinition->get(name => 'fastqs_split', type => 'fastq', output_sub => sub { my ($self, $step) = @_;
                                                                                                                         my $input_files = $step->inputs->{fastq_files};
                                                                                                                         my @fastq_files;
                                                                                                                         foreach my $vrfile (@$input_files) {
                                                                                                                             push(@fastq_files, $vrfile->path->basename.'.md5');
                                                                                                                         }
                                                                                                                         return [@fastq_files]; }) };
    }
    method description {
        return "Takes a single-ended fastq file and/or 2 fastq files of a paired-end read and splits them into multiple smaller fastq files";
    }
}

1;
