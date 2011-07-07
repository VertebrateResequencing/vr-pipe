use VRPipe::Base;

role VRPipe::StepRole {
    method name {
        my $class = ref($self);
        my ($name) = $class =~ /.+::(.+)/;
        return $name;
    }
    requires 'options_definition';
    requires 'inputs_definition';
    requires 'body_sub';
    requires 'post_process_sub';
    requires 'outputs_definition';
    requires 'description';
    
    # these may be needed by body_sub and post_process_sub
    has 'step_state' => (is => 'rw',
                         isa => 'VRPipe::StepState');
 
    has 'data_element' => (is => 'ro',
                           isa => 'VRPipe::DataElement',
                           builder => '_build_data_element',
                           lazy => 1);
    
    has 'output_root' => (is => 'ro',
                          isa => Dir,
                          coerce => 1,
                          builder => '_build_output_root',
                          lazy => 1);
    
    has 'options' => (is => 'ro',
                      isa => 'HashRef',
                      builder => '_resolve_options',
                      lazy => 1);
    
    has 'inputs' => (is => 'ro',
                     isa => PersistentFileHashRef,
                     builder => '_resolve_inputs',
                     lazy => 1);
    
    has 'outputs' => (is => 'ro',
                      isa => PersistentFileHashRef,
                      builder => '_build_outputs',
                      lazy => 1);
    has 'temps' => (is => 'ro',
                    isa => ArrayRefOfPersistent,
                    builder => '_build_temps',
                    lazy => 1);
    
    has 'previous_step_outputs' => (is => 'rw',
                                    isa => PersistentFileHashRef);
    
    # when parse is called, we'll store our dispatched refs here
    has 'dispatched' => (is => 'ro',
                         traits  => ['Array'],
                         isa     => 'ArrayRef', #*** ArrayRef[ArrayRef[Str,VRPipe::Requirements]] doesn't work, don't know why...
                         lazy    => 1,
                         default => sub { [] },
                         handles => { dispatch => 'push',
                                      num_dispatched  => 'count' });
    
    # and we'll also store all the output files the body_sub makes
    has '_output_files' => (is => 'ro',
                            traits  => ['Hash'],
                            isa     => 'HashRef',
                            lazy    => 1,
                            default => sub { {} },
                            handles => { _remember_output_files => 'set' });
    has '_temp_files' => (is => 'ro',
                          traits  => ['Array'],
                          isa     => 'ArrayRef',
                          lazy    => 1,
                          default => sub { [] },
                          handles => { _remember_temp_file => 'push' });
    
    method _build_data_element {
        my $step_state = $self->step_state || $self->throw("Cannot get data element without step state");
        return $step_state->dataelement;
    }
    method _build_output_root {
        my $step_state = $self->step_state || $self->throw("Cannot get output root without step state");
        return $step_state->pipelinesetup->output_root;
    }
    method _build_outputs {
        my $step_state = $self->step_state || $self->throw("Cannot get outputs without step state");
        return $step_state->output_files;
    }
    method _build_temps {
        my $step_state = $self->step_state || $self->throw("Cannot get outputs without step state");
        return $step_state->temp_files;
    }
    
    method _resolve_options {
        my $step_state = $self->step_state || $self->throw("Cannot get options without step state");
        my $user_opts = $step_state->pipelinesetup->options;
        my $hash = $self->options_definition;
        
        my %return;
        while (my ($key, $val) = each %$hash) {
            if ($val->isa('VRPipe::StepOption')) {
                my $user_val = $user_opts->{$key};
                if (defined $user_val) {
                    my $allowed = $val->allowed_values;
                    if (@$allowed) {
                        my %allowed = map { $_ => 1 } @$allowed;
                        if (exists $allowed{$user_val}) {
                            $return{$key} = $user_val;
                        }
                        else {
                            $self->throw("'$user_val' is not an allowed option for '$key'");
                        }
                    }
                    else {
                        $return{$key} = $user_val;
                    }
                }
                elsif (! $val->optional) {
                    $self->throw("the option '$key' is required");
                }
                else {
                    my $default = $val->default_value;
                    if (defined $default) {
                        $return{$key} = $default;
                    }
                }
            }
            else {
                $self->throw("invalid class ".ref($val)." supplied for option '$key' definition");
            }
        }
        
        return \%return;
    }
    
    method _resolve_inputs {
        my $hash = $self->inputs_definition;
        my $step_num = $self->step_state->stepmember->step_number;
        my $step_adaptor = VRPipe::StepAdaptor->get(pipeline => $self->step_state->stepmember->pipeline, before_step_number => $step_num);
        
        my %return;
        while (my ($key, $val) = each %$hash) {
            if ($val->isa('VRPipe::File')) {
                $return{$key} = [$val];
            }
            elsif ($val->isa('VRPipe::StepIODefinition')) {
                # see if we have this $key in our previous_step_outputs or
                # via the options or data_element
                my $results;
                
                my $pso = $self->previous_step_outputs;
                if ($pso) {
                    # can our StepAdaptor adapt previous step output to this
                    # input key?
                    $results = $step_adaptor->adapt(input_key => $key, pso => $pso);
                }
                if (! $results) {
                    # can our StepAdaptor translate our dataelement into a file
                    # for this key?
                    $results = $step_adaptor->adapt(input_key => $key, data_element => $self->data_element);
                }
                if (! $results) {
                    my $opts = $self->options;
                    if ($opts && defined $opts->{$key}) {
                        #*** the step should have had some way of advertising
                        #    what its inputs_definition keys were, so that the
                        #    user could have provided values for them during
                        #    PipelineSetup
                        $results = [VRPipe::File->get(path => $opts->{$key})];
                    }
                }
                
                if (! $results) {
                    $self->throw("the input file(s) for '$key' of stepstate ".$self->step_state->id." could not be resolved");
                }
                
                my $num_results = @$results;
                my $max_allowed = $val->max_files;
                my $min_allowed = $val->min_files;
                if ($max_allowed == -1) {
                    $max_allowed = $num_results;
                }
                if ($min_allowed == -1) {
                    $min_allowed = $num_results;
                }
                unless ($num_results >= $min_allowed && $num_results <= $max_allowed) {
                    $self->throw("there were $num_results input file(s) for '$key' of stepstate ".$self->step_state->id.", which does not fit the allowed range $min_allowed..$max_allowed");
                }
                
                my @vrfiles;
                foreach my $result (@$results) {
                    unless (ref($result) && ref($result) eq 'VRPipe::File') {
                        $result = VRPipe::File->get(path => file($result)->absolute);
                    }
                    
                    my $type = VRPipe::FileType->create($val->type, {file => $result->path});
                    unless ($type->check_type) {
                        $self->throw("file ".$result->path." was not the correct type");
                    }
                    
                    push(@vrfiles, $result);
                }
                
                $return{$key} = \@vrfiles;
            }
            else {
                $self->throw("invalid class ".ref($val)." supplied for input '$key' value definition");
            }
        }
        
        return \%return;
    }
    
    method _missing (PersistentFileHashRef $hash, PersistentHashRef $defs) {
        # check that we don't have any outputs defined in the definition that
        # no files were made for
        while (my ($key, $val) = each %$defs) {
            next if exists $hash->{$key};
            $self->throw("'$key' was defined as an output, yet no output file was made with that output_key");
        }
        
        my @missing;
        # check the files we actually output are as expected
        while (my ($key, $val) = each %$hash) {
            foreach my $file (@$val) {
                if (! $file->s) {
                    push(@missing, $file->path);
                }
                else {
                    my $bad = 0;
                    
                    # check the filetype is correct
                    my $type = VRPipe::FileType->create($file->type, {file => $file->path});
                    unless ($type->check_type) {
                        $self->warn($file->path." exists, but is the wrong type!");
                        $bad = 1;
                    }
                    
                    # check the expected metadata keys exist
                    my $def = $defs->{$key};
                    if ($def && $def->isa('VRPipe::StepIODefinition')) {
                        my @needed = $def->required_metadata_keys;
                        if (@needed) {
                            my $meta = $file->metadata;
                            foreach my $key (@needed) {
                                unless (exists $meta->{$key}) {
                                    $self->warn($file->path." exists, but lacks metadata key $key!");
                                    $bad = 1;
                                }
                            }
                        }
                    }
                    
                    if ($bad) {
                        push(@missing, $file->path);
                    }
                }
            }
        }
        
        return @missing;
    }
    
    method missing_input_files {
        return $self->_missing($self->inputs, $self->inputs_definition);
    }
    
    method missing_output_files {
        return $self->_missing($self->outputs, $self->outputs_definition);
    }
    
    method _run_coderef (Str $method_name) {
        my $ref = $self->$method_name();
        return &$ref($self);
    }
    
    method output_file (Str :$output_key, File|Str :$basename, FileType :$type, Dir|Str :$output_dir?, HashRef :$metadata?, Bool :$temporary = 0) {
        $output_dir ||= $self->output_root;
        
        my $vrfile = VRPipe::File->get(path => file($output_dir, $basename), type => $type);
        $vrfile->add_metadata($metadata) if $metadata;
        
        my $hash = $self->_output_files;
        my $files = $hash->{$output_key} || [];
        push(@$files, $vrfile);
        $self->_remember_output_files($output_key => $files);
        
        if ($temporary) {
            $self->_remember_temp_file($vrfile);
        }
        
        return $vrfile;
    }
    
    method parse {
        my @missing = $self->missing_input_files;
        $self->throw("Required input files are missing: (@missing)") if @missing;
        
        $self->_run_coderef('body_sub');
        
        # store output and temp files on the StepState
        my $output_files = $self->_output_files;
        if (keys %$output_files) {
            my $step_state = $self->step_state;
            $step_state->output_files($output_files);
            $step_state->update;
        }
        my $temp_files = $self->_temp_files;
        if (@$temp_files) {
            my $step_state = $self->step_state;
            $step_state->temp_files($temp_files);
            $step_state->update;
        }
        
        # return true if we've finished the step
        my $dispatched = $self->num_dispatched;
        if ($dispatched) {
            return 0;
        }
        else {
            return $self->post_process;
        }
    }
    
    method post_process {
        my $ok = $self->_run_coderef('post_process_sub');
        my $debug_desc = "step ".$self->name." failed for data element ".$self->data_element->id." and pipelinesetup ".$self->step_state->pipelinesetup->id;
        my $stepstate = $self->step_state;
        if ($ok) {
            my @missing = $self->missing_output_files;
            $stepstate->unlink_temp_files;
            if (@missing) {
                $self->throw("Some output files are missing (@missing) for $debug_desc");
            }
            else {
                return 1;
            }
        }
        else {
            $stepstate->unlink_temp_files;
            $self->throw("The post-processing part of $debug_desc");
        }
    }
    
    method new_requirements (Int :$memory, Int :$time, Int :$cpus?, Int :$tmp_space?, Int :$local_space?, HashRef :$custom?) {
        return VRPipe::Requirements->get(memory => $memory,
                                         time => $time,
                                         $cpus ? (cpus => $cpus) : (),
                                         $tmp_space ? (tmp_space => $tmp_space) : (),
                                         $local_space ? (local_space => $local_space) : (),
                                         $custom ? (custom => $custom) : ());
    }
    
    method dispatch_vrpipecode (Str $code, VRPipe::Requirements $req) {
        my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
        my $cmd = qq[perl -MVRPipe::Persistent::Schema -e "VRPipe::Persistent::SchemaBase->database_deployment(q[$deployment]); $code"];
        $self->dispatch([$cmd, $req]);
    }
    
    method dispatch_md5sum (VRPipe::File $vrfile, Maybe[Str] $expected_md5) {
        my $path = $vrfile->path;
        my $req = $self->new_requirements(memory => 50, time => 1);
        
        if ($expected_md5) {
            return $self->dispatch_vrpipecode(qq[use Digest::MD5; open(FILE, q[$path]) or die q[Could not open file $path]; binmode(FILE); if (Digest::MD5->new->addfile(*FILE)->hexdigest eq q[$expected_md5]) { VRPipe::File->get(path => q[$path], md5 => q[$expected_md5]); } else { die q[md5sum of $path does not match expected value] }],
                                              $req);
        }
        else {
            return $self->dispatch_vrpipecode(qq[use Digest::MD5; open(FILE, q[$path]) or die q[Could not open file $path]; binmode(FILE); VRPipe::File->get(path => q[$path], md5 => Digest::MD5->new->addfile(*FILE)->hexdigest);],
                                              $req);
        }
    }
}

1;