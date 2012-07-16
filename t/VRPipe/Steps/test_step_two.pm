use VRPipe::Base;

class VRPipe::Steps::test_step_two with VRPipe::StepRole {
    method options_definition {
        return { all_option => VRPipe::StepOption->create(description => 'an option that applies to all steps'),
                 two_option => VRPipe::StepOption->create(description => 'an optional option for step two',
                                                       optional => 1) };
    }
    method inputs_definition {
        return { two_input => VRPipe::StepIODefinition->create(type => 'txt',
                                                            description => 'step two input file',
                                                            metadata => {one_meta => 'metadata we require to appear on our input file'}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $all_opt = Path::Class::File->new($options->{all_option});
            my $two_opt = $options->{two_option} || 'body_decided_two_option';
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'cat',
                                                               version => VRPipe::StepCmdSummary->determine_version('cat --version', '^cat \(GNU coreutils\) (\S+)$'),
                                                               summary => 'cat $input_file > $output_file'));
            
            my $req = $self->new_requirements(memory => 50, time => 1);
            foreach my $in (@{$self->inputs->{two_input}}) {
                my $out = $self->output_file(output_key => 'two_output',
                                             basename => $in->basename.'.step_two',
                                             type => 'txt',
                                             metadata => {one_meta => $in->metadata->{one_meta},
                                                          two_meta => $two_opt});
                
                my ($in_path, $out_path) = ($in->path, $out->path);
                my $this_cmd = "cat $in_path > $out_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::test_step_two', 'do_step_two', [$this_cmd, $req]);
            }
        };
    }
    method outputs_definition {
        return { two_output => VRPipe::StepIODefinition->create(type => 'txt',
                                                             description => 'step two output file',
                                                             metadata => {one_meta => 'metadata from step one',
                                                                          two_meta => 'metadata applied to step two output file'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "test step two";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method do_step_two (ClassName|Object $self: Str $cmd_line) {
        # do some pre-processing...
        system($cmd_line) && $self->throw("[$cmd_line] failed");
        # do some post-processing...
    }
}

1;