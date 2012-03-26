use VRPipe::Base;

class VRPipe::Steps::java with VRPipe::StepRole {
    use POSIX qw(ceil);
    
    has 'java_exe' => (is => 'rw',
                       isa => 'Str');
    
    has 'memory_multiplier' => (is => 'ro',
                                isa => 'Num',
                                default => 0.9);
    
    has 'standard_options' => (is => 'ro',
                               isa => 'ArrayRef',
                               lazy => 1,
                               builder => '_build_standard_options');
    
    method _build_smaller_recommended_requirements_override {
        return 0;
    }
    
    method _build_standard_options {
        return ['java_exe'];
    }
    
    method jvm_args (ClassName|Object $self: Int $memory, Str|Dir $dir?) {
        my $java_mem = ceil($memory * $self->memory_multiplier);
        if ($java_mem < 50) {
            $java_mem = 50;
        }
        
        my $temp_dir = '';
        if ($dir) {
            $temp_dir = ' -Djava.io.tmpdir='.$self->tempdir(DIR => $dir);
        }
        return qq[-Xmx${java_mem}m -Xms${java_mem}m$temp_dir -server -XX:+UseSerialGC];
    }
    
    method handle_standard_options (HashRef $options) {
        foreach my $method (@{$self->standard_options}) {
            next unless defined $options->{$method};
            $self->$method($options->{$method});
        }
    }
    
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
