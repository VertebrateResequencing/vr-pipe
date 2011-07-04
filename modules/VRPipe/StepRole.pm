use VRPipe::Base;

role VRPipe::StepRole {
    method name {
        my $class = ref($self);
        my ($name) = $class =~ /.+::(.+)/;
        return $name;
    }
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
                      builder => '_build_options',
                      lazy => 1);
    
    has 'inputs' => (is => 'ro',
                     isa => PersistentFileHashRef,
                     builder => 'resolve_inputs',
                     lazy => 1);
    
    has 'outputs' => (is => 'ro',
                      isa => PersistentFileHashRef,
                      builder => '_build_outputs',
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
    
    method _build_data_element {
        my $step_state = $self->step_state || $self->throw("Cannot get data element without step state");
        return $step_state->dataelement;
    }
    method _build_output_root {
        my $step_state = $self->step_state || $self->throw("Cannot get output root without step state");
        return $step_state->pipelinesetup->output_root;
    }
    method _build_options {
        my $step_state = $self->step_state || $self->throw("Cannot get options without step state");
        return $step_state->pipelinesetup->options;
    }
    method _build_outputs {
        my $step_state = $self->step_state || $self->throw("Cannot get outputs without step state");
        return $step_state->output_files;
    }
    
    method resolve_inputs {
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
        my @missing;
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
                    if ($def->isa('VRPipe::StepIODefinition')) {
                        my $def_meta = $def->metadata;
                        my @needed = keys %$def_meta;
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
    
    method output_file (Str :$output_key, File|Str :$basename, FileType :$type, Dir|Str :$output_dir?, HashRef :$metadata?) {
        $output_dir ||= $self->output_root;
        
        my $vrfile = VRPipe::File->get(path => file($output_dir, $basename), type => $type);
        $vrfile->add_metadata($metadata) if $metadata;
        
        my $hash = $self->_output_files;
        my $files = $hash->{$output_key} || [];
        push(@$files, $vrfile);
        $self->_remember_output_files($output_key => $files);
        
        return $vrfile;
    }
    
    method parse {
        my @missing = $self->missing_input_files;
        $self->throw("Required input files are missing: (@missing)") if @missing;
        
        $self->_run_coderef('body_sub');
        
        # store output files on the StepState
        my $output_files = $self->_output_files;
        if (keys %$output_files) {
            my $step_state = $self->step_state;
            $step_state->output_files($output_files);
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
        # in case our main coderef created a file without using
        # VRPipe::File->openw and ->close, we have to force a stats update for
        # all output files
        my $outputs = $self->outputs;
        if ($outputs) {
            foreach my $val (values %$outputs) {
                foreach my $file (@$val) {
                    if (! $file->e || ! $file->s) {
                        $file->update_stats_from_disc;
                    }
                }
            }
        }
        
        my $ok = $self->_run_coderef('post_process_sub');
        my $debug_desc = "step ".$self->name." failed for data element ".$self->data_element->id." and pipelinesetup ".$self->step_state->pipelinesetup->id;
        if ($ok) {
            my @missing = $self->missing_output_files;
            if (@missing) {
                $self->throw("Some output files are missing (@missing) for $debug_desc");
            }
            else {
                #*** concept of output paths marked as temporary; delete them
                #    now
                
                return 1;
            }
        }
        else {
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
}

1;