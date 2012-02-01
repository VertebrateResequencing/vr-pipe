use VRPipe::Base;

class VRPipe::Steps::cramtools extends VRPipe::Steps::java {
    has 'cramtools_path' => (is => 'rw',
                             isa => Dir,
                             coerce => 1);
    
    around _build_standard_options {
        return [@{$self->$orig}, 'cramtools_path'];
    }
    
    our %CRAMTOOLS_VERSIONS;
    has 'cramtools_version' => (is => 'ro',
                                isa => 'Str',
                                lazy => 1,
                                builder => 'determine_cramtools_version');
    
    method jar (ClassName|Object $self:) {
        return file($self->cramtools_path, 'cramtools.jar');
    }
    
    method determine_cramtools_version (ClassName|Object $self:) {
        my $cramtools_jar = $self->jar->stringify;
        unless (defined $CRAMTOOLS_VERSIONS{$cramtools_jar}) {
            my $jvm_args = $self->jvm_args(50);
            my $java_exe = $self->java_exe;
            $CRAMTOOLS_VERSIONS{$cramtools_jar} = VRPipe::StepCmdSummary->determine_version(qq[$java_exe $jvm_args -jar $cramtools_jar -h], 'v([\d\.\-]+[a-z\d]*)');
        }
        return $CRAMTOOLS_VERSIONS{$cramtools_jar};
    }
    
    around options_definition {
        return { %{$self->$orig},
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
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
        return "Generic step for using the cramtools";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
