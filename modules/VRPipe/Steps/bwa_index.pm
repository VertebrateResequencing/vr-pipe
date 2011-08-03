use VRPipe::Base;

class VRPipe::Steps::bwa_index with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file to map against'),
                 bwa_index_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa index command line, including desired options, but excluding the reference fasta file',
                                                          optional => 1,
                                                          default_value => 'bwa index -a bwtsw') };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $ref = Path::Class::File->new($self->options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $cmd = $self->options->{bwa_index_cmd};
            my @c = split(" ", $cmd);
            unless ($c[1] eq 'index') {
                $self->throw("bad bwa_index_cmd '$cmd'");
            }
            
            if ($cmd =~ /$ref/) {
                $self->throw("bwa_index_cmd should not include the reference");
            }
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => Path::Class::File->new($c[0])->basename, version => VRPipe::StepCmdSummary->determine_version($c[0], '^Version: (.+)$'), summary => $cmd.' $reference_fasta'));
            $cmd .= ' '.$ref;
            
            foreach my $suffix (qw(bwt pac rbwt rpac rsa sa)) {
                $self->output_file(output_key => 'bwa_index_binary_files', output_dir => $ref->dir->stringify, basename => $ref->basename.'.'.$suffix, type => 'bin');
            }
            foreach my $suffix (qw(amb ann)) {
                $self->output_file(output_key => 'bwa_index_text_files', output_dir => $ref->dir->stringify, basename => $ref->basename.'.'.$suffix, type => 'txt');
            }
            
            $self->dispatch([$cmd, $self->new_requirements(memory => 3900, time => 1), {block_and_skip_if_ok => 1}]);
        };
    }
    method outputs_definition {
        return { bwa_index_binary_files => VRPipe::StepIODefinition->get(type => 'bin', description => 'the files produced by bwa index', min_files => 6, max_files => 6),
                 bwa_index_text_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'the files produced by bwa index', min_files => 2, max_files => 2) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Indexes a reference genome fasta file, making it suitable for use in subsequent bwa mapping";
    }
}

1;
