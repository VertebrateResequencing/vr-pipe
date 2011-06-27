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
                      builder => 'resolve_outputs',
                      lazy => 1);
    
    has 'previous_step_outputs' => (is => 'rw',
                                    isa => PersistentFileHashRef);
    
    # when parse is called, we'll store our dispatched refs here
    has 'dispatched' => (is => 'ro',
                         traits  => ['Array'],
                         isa     => 'ArrayRef', #*** ArrayRef[ArrayRef[Str,VRPipe::Requirements]] doesn't work, don't know why...
                         lazy    => 1,
                         default => sub { [] },
                         handles => { dispatch => 'push' });
    
    ## prior to derefrencing and running the step subroutine, run_step does:
    #my @missing = $obj->missing_required_files($step_name);
    ## if any files are missing, then run_step returns false and does nothing.
    ## Otherwise it runs the desired subroutine. If the subroutine returns true,
    ## that means the step completed and we can move on to the finish_step
    ## step. If it returned false, it probably used dispatch() (see below) so we
    ## see what it dispatched:
    #my @dispatches = $obj->dispatched($step_name);
    ## this returns a list of [$cmd_line_string, $requirements_object] refs that
    ## you could use to create VRPipe::JobManager::Submission objects and
    ## eventually actually get the cmd_line to run. You'll keep records on how
    ## those Submissions do, and when they're all successfull you'll do:
    #$obj->finish_step($step_name);
    ## this runs a post-processing method for the step (that is not allowed
    ## to do any more dispatching) and returns true if the post-process was fine.
    ## finish_step appends its own auto-check that all the provided files
    ## exist. So if finish_step returns true you record in your system that
    ## this step (on this pipeline setup and datasource element) finished, so
    ## the next time you come to this loop you'll skip and try the following
    ## step. (VRPipe::PipelineManager::Manager does all this sort of stuff for
    ## you.)
    
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
        return $step_state->pipelinesetup->options; #*** supposed to be a hash ref, currently a str...
    }
    
    method resolve_inputs {
        my $hash = $self->inputs_definition;
        
        my %return;
        while (my ($key, $val) = each %$hash) {
            if ($val->isa('VRPipe::File')) {
                $return{$key} = $val;
            }
            elsif ($val->isa('VRPipe::FileDefinition')) {
                # see if we have this $key in our previous_step_outputs or
                # via the options or data_element
                my $input_vrfile;
                my @vrpfile_get_args = $val->type eq 'any' ? () : (type => $val->type); #*** need a proper way of not setting generic types if file already has a specific one
                my $pso = $self->previous_step_outputs;
                if ($pso) {
                    #*** StepAdaptor for connecting output keys of one step with
                    #    input keys of another step, without having to alter the
                    #    io definitions of either step...
                    if (defined $pso->{$key} && $val->matches($pso->{$key})) {
                        $input_vrfile = $pso->{$key};
                    }
                }
                if (! $input_vrfile) {
                    my $der = $self->data_element->result;
                    #*** StepAdaptor adapts data_elements as well?
                    #    Don't yet know how to actually get a file out of a DE,
                    #    if the DE even represents files...
                    if ($der && $val->matches($der)) {
                        $input_vrfile = VRPipe::File->get(path => file($der)->absolute, @vrpfile_get_args);
                    }
                }
                if (! $input_vrfile) {
                    my $opts = $self->options;
                    if ($opts) {
                        #*** the step should have had some way of advertising
                        #    what its inputs_definition keys were, so that the
                        #    user could have provided values for them during
                        #    PipelineSetup
                        if (defined $opts->{$key} && $val->matches($opts->{$key})) {
                            $input_vrfile = VRPipe::File->get(path => $opts->{$key}, @vrpfile_get_args);
                        }
                    }
                }
                
                if (! $input_vrfile) {
                    $self->throw("the input file for '$key' of stepstate ".$self->step_state->id." could not be resolved");
                }
                
                $return{$key} = $input_vrfile;
            }
            else {
                $self->throw("invalid class ".ref($val)." supplied for input '$key' value definition");
            }
        }
        
        return \%return;
    }
    
    method resolve_outputs {
        my $hash = $self->outputs_definition;
        my $output_root = $self->output_root;
        
        my %return;
        while (my ($key, $val) = each %$hash) {
            if ($val->isa('VRPipe::File')) {
                $return{$key} = $val;
            }
            elsif ($val->isa('VRPipe::FileDefinition')) {
                my $basename = $val->output_basename($self);
                if (ref($basename) eq 'ARRAY') {
                    my @vrf_array;
                    foreach my $bn (@$basename) {
                        push(@vrf_array, VRPipe::File->get(path => file($output_root, $bn), $val->type eq 'any' ? () : (type => $val->type)));
                    }
                    $return{$key} = \@vrf_array;
                }
                else {
                    $return{$key} = VRPipe::File->get(path => file($output_root, $basename), type => $val->type);
                }
            }
            else {
                $self->throw("invalid class ".ref($val)." supplied for output '$key' value definition");
            }
        }
        
        return \%return;
    }
    
    method _missing (PersistentFileHashRef $hash, Bool $check_type) {
        my @missing;
        while (my ($key, $val) = each %$hash) {
            my @files;
            if (ref($val) eq 'ARRAY') {
                @files = @$val;
            }
            else {
                @files = ($val);
            }
            
            foreach my $file (@files) {
                if (! $file->s) {
                    push(@missing, $file->path);
                }
                elsif ($check_type) {
                    my $type = VRPipe::FileType->create($file->type, {file => $file->path});
                    unless ($type->check_type) {
                        $self->warn($file->path." exists, but is the wrong type!");
                        push(@missing, $file->path);
                    }
                }
            }
        }
        return @missing;
    }
    
    method missing_input_files {
        return $self->_missing($self->inputs, 0);
    }
    
    method missing_output_files {
        return $self->_missing($self->outputs, 1);
    }
    
    method _run_coderef (Str $method_name) {
        my $ref = $self->$method_name();
        return &$ref($self);
    }
    
    method parse {
        my @missing = $self->missing_input_files;
        $self->throw("Required input files are missing: (@missing)") if @missing;
        
        my $finished = $self->_run_coderef('body_sub');
        if ($finished) {
            return $self->post_process;
        }
        
        # presumably body_sub called dispatch(), so we'll leave it up to our
        # caller to do stuff with dispatched()
        return 0;
    }
    
    method post_process {
        # in case our main coderef created a file without using
        # VRPipe::File->openw and ->close, we have to force a stats update for
        # all output files
        my $outputs = $self->outputs;
        if ($outputs) {
            foreach my $val (values %$outputs) {
                my @files;
                if (ref($val) eq 'ARRAY') {
                    @files = @$val;
                }
                else {
                    @files = ($val);
                }
                
                foreach my $file (@files) {
                    if (! $file->e || ! $file->s) {
                        $file->update_stats_from_disc; #*** when can we force an md5 check? Not here since outside of bsub...
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
    
    method _input_files_from_input_or_element (Str $input_key) {
        my $input = $self->inputs->{$input_key} || $self->data_element->result;
        my @inputs;
        if (ref($input) eq 'ARRAY') {
            @inputs = @{$input};
        }
        else {
            @inputs = ($input);
        }
        my @files;
        foreach my $input (@inputs) {
            if (ref($input) && $input->isa('VRPipe::File')) {
                push(@files, $input->path);
            }
            else {
                push(@files, file($input));
            }
        }
        return [@files];
    }
}

1;