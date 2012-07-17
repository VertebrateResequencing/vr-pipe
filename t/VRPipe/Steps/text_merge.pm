use VRPipe::Base;

class VRPipe::Steps::text_merge with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { input_text_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1,
                                                                   description => 'text files to merge') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'cat',
                                                               version => VRPipe::StepCmdSummary->determine_version('cat --version', '^cat \(GNU coreutils\) (\S+)$'),
                                                               summary => 'cat $input_file(s) > $output_file'));
            
            my $req = $self->new_requirements(memory => 50, time => 1);
            my @in_files = @{$self->inputs->{input_text_files}};
            my @in_paths = map { $_->path} @in_files;
            my $out = $self->output_file(output_key => 'merged_file',
                                         basename => 'merged.txt',
                                         type => 'txt',
                                         metadata => { source_files => join ',', @in_paths });
                
            my $this_cmd = "cat @in_paths > ".$out->path;
                
            $self->dispatch([$this_cmd, $req, {output_files => [$out]}]);
        };
    }
    method outputs_definition {
        return { merged_file => VRPipe::StepIODefinition->create(type => 'txt',
                                                              description => 'merged text file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "merges text files into a single file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;