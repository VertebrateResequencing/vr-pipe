use VRPipe::Base;

class VRPipe::Steps::java with VRPipe::StepRole {
    method options_definition {
        return { java_exe => VRPipe::StepOption->get(description => 'path to your java executable', optional => 1, default_value => 'java'),
                 tmp_dir => VRPipe::StepOption->get(description => 'location for tmp directories; defaults to working directory', optional => 1),
                };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub { return 1; };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generic step for steps using java";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
