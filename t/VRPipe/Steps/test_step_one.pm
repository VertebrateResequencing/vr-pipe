use VRPipe::Base;

class VRPipe::Steps::test_step_one with VRPipe::StepRole {
    method options_definition {
        return { all_option => VRPipe::StepOption->create(description => 'an option that applies to all steps'),
                 one_option => VRPipe::StepOption->create(description => 'a required option for step one') };
    }
    
    method inputs_definition {
        return { one_input => VRPipe::StepIODefinition->create(type        => 'txt',
                                                               description => 'step one input file') };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $all_opt = Path::Class::File->new($options->{all_option});
            my $one_opt = $options->{one_option};
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe     => 'cat',
                                                                  version => VRPipe::StepCmdSummary->determine_version('cat --version', '^cat \(GNU coreutils\) (\S+)$'),
                                                                  summary => 'cat $input_file > $output_file'));
            
            my $req = $self->new_requirements(memory => 50, time => 1);
            foreach my $in (@{ $self->inputs->{one_input} }) {
                my $out = $self->output_file(output_key => 'one_output',
                                             basename   => $in->basename . '.step_one',
                                             type       => 'txt',
                                             metadata   => { one_meta => $one_opt });
                
                my ($in_path, $out_path) = ($in->path, $out->path);
                my $this_cmd = "cat $in_path > $out_path";
                
                $self->dispatch([$this_cmd, $req]);
            }
        };
    }
    
    method outputs_definition {
        return { one_output => VRPipe::StepIODefinition->create(type        => 'txt',
                                                                description => 'step one output file',
                                                                metadata    => { one_meta => 'metadata applied to step one output file' }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "test step one";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
