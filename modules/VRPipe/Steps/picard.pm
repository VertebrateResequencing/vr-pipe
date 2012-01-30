use VRPipe::Base;

class VRPipe::Steps::picard extends VRPipe::Steps::java {
    has 'picard_path' => (is => 'rw',
                          isa => Dir,
                          coerce => 1);
    
    has '+memory_multiplier' => (default => 0.7);
    
    around _build_standard_options {
        return [@{$self->$orig}, 'picard_path'];
    }
    
    our %PICARD_VERSIONS;
    has 'picard_version' => (is => 'ro',
                             isa => 'Str',
                             lazy => 1,
                             builder => 'determine_picard_version');
    
    method determine_picard_version (ClassName|Object $self:) {
        my $picard_path = $self->picard_path->stringify;
        unless (defined $PICARD_VERSIONS{$picard_path}) {
            my $version = 0;
            opendir(my $dh, $picard_path) || $self->throw("Could not open picard directory $picard_path");
            foreach (readdir $dh) 
            {
                if (/^picard-([\d\.]+)\.jar/)
                {
                    $version = $1;
                    last;
                }
            }
            closedir($dh);
            $PICARD_VERSIONS{$picard_path} = $version;
        }
        return $PICARD_VERSIONS{$picard_path};
    }
    
    method jar (ClassName|Object $self: Str $basename!) {
        return file($self->picard_path, $basename);
    }
    
    around options_definition {
        return { %{$self->$orig},
                 picard_path => VRPipe::StepOption->get(description => 'path to Picard jar files', optional => 1, default_value => "$ENV{PICARD}"),
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
        return "Generic step for using Picard";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
