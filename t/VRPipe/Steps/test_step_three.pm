use VRPipe::Base;

class VRPipe::Steps::test_step_three with VRPipe::StepRole {
    method options_definition {
        return { all_option   => VRPipe::StepOption->create(description => 'an option that applies to all steps'),
                 three_option => VRPipe::StepOption->create(description   => 'an optional option for step three with a default',
                                                            optional      => 1,
                                                            default_value => 'StepOption_default_decided_three_option') };
    }
    
    method inputs_definition {
        return { three_input => VRPipe::StepIODefinition->create(type        => 'txt',
                                                                 description => 'step three input file',
                                                                 metadata    => {
                                                                               one_meta => 'metadata we require to appear on our input file, from step one',
                                                                               two_meta => 'metadata we require to appear on our input file, from step two' }) };
    }
    
    method body_sub {
        return sub {
            my $self      = shift;
            my $options   = $self->options;
            my $all_opt   = Path::Class::File->new($options->{all_option});
            my $three_opt = $options->{three_option};
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe     => 'cat',
                                                                  version => VRPipe::StepCmdSummary->determine_version('cat --version', '^cat \(GNU coreutils\) (\S+)$'),
                                                                  summary => 'cat $input_file > $output_file'));
            
            my $req = $self->new_requirements(memory => 50, time => 1);
            foreach my $in (@{ $self->inputs->{three_input} }) {
                my $out = $self->output_file(output_key => 'three_output',
                                             basename   => $in->basename . '.step_three',
                                             type       => 'txt',
                                             metadata   => {
                                                           one_meta => $in->metadata->{one_meta},
                                                           two_meta => $in->metadata->{two_meta},
                                                           $in->basename =~ /^file2/ ? (three_meta => $three_opt) : () });
                
                my ($in_path, $out_path) = ($in->path, $out->path);
                
                $self->dispatch_vrpipecode(qq[use VRPipe::Steps::test_step_three; VRPipe::Steps::test_step_three->do_step_three(q[$in_path], q[$out_path])], $req);
            }
        };
    }
    
    method outputs_definition {
        return { three_output => VRPipe::StepIODefinition->create(type        => 'txt',
                                                                  description => 'step three output file',
                                                                  metadata    => {
                                                                                one_meta   => 'metadata from step one',
                                                                                two_meta   => 'metadata from step two',
                                                                                three_meta => 'metadata that might be applied to step three output file',
                                                                                optional   => ['three_meta'] }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "test step three";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method do_step_three (ClassName|Object $self: Str|File $in_path, Str|File $out_path) {
        # do some complicated custom perl, not possible with a simple cmd line...
        my $cmd_line = "cat $in_path > $out_path";
        system($cmd_line) && $self->throw("[$cmd_line] failed");
    }
}

1;
