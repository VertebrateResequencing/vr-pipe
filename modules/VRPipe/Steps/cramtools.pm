use VRPipe::Base;

class VRPipe::Steps::cramtools extends VRPipe::Steps::java {
    around options_definition {
        return { %{$self->$orig},
                 cramtools_path => VRPipe::StepOption->get(description => 'path to cramtools jar file', optional => 1, default_value => "$ENV{CRAMTOOLS}"),
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
        return "Generic step for using cramtools";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
