use VRPipe::Base;

class VRPipe::Step extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'name' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    # *** the following attributes are supposed to store code refs...
    
    has 'inputs_sub' => (is => 'rw',
                         isa => 'CodeRef',
                         traits => ['VRPipe::Persistent::Attributes']);
    
    has 'body_sub' => (is => 'rw',
                       isa => 'CodeRef',
                       traits => ['VRPipe::Persistent::Attributes']);
    
    has 'post_process_sub' => (is => 'rw',
                               isa => 'CodeRef',
                               traits => ['VRPipe::Persistent::Attributes']);
    
    has 'outputs_sub' => (is => 'rw',
                          isa => 'CodeRef',
                          traits => ['VRPipe::Persistent::Attributes']);
    
    has 'description' => (is => 'rw',
                         isa => Varchar[64],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_nullable => 1);
    
    # these transients may be needed by inputs, body, post_process and outputs
    # code refs to do their jobs
    has 'step_state' => (is => 'rw',
                         isa => 'VRPipe::StepState');
 
    has 'data_element' => (is => 'ro',
                           isa => 'VRPipe::DataElement',
                           builder => '_build_data_element',
                           lazy => 1);
    
    has 'output_root' => (is => 'ro',
                          isa => Dir,
                          builder => '_build_output_root',
                          lazy => 1);
    
    has 'options' => (is => 'ro',
                      isa => 'HashRef',
                      builder => '_build_options',
                      lazy => 1);
    
    has 'input' => (is => 'rw',
                    isa => 'Defined');
    
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
    
    __PACKAGE__->make_persistent();
    
    
    method _run_coderef (Str $method_name) {
        my $ref = $self->$method_name();
        return &$ref($self);
    }
    
    method missing_input_files {
        my @missing;
        foreach my $input ($self->_run_coderef('inputs_sub')) {
            #*** should use db for checking file existance
            if (! -s $input) {
                push(@missing, $input);
            }
        }
        
        return @missing;
    }
    
    method missing_output_files {
        my @missing;
        foreach my $output ($self->_run_coderef('outputs_sub')) {
            #*** should use db for checking file existance
            if (! -s $output) {
                push(@missing, $output);
            }
        }
        
        return @missing;
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
}

1;