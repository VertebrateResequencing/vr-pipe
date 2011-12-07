use VRPipe::Base;

role VRPipe::PipelineRole {
    use VRPipe::StepAdaptorDefiner;
    use VRPipe::StepBehaviourDefiner;
    use Module::Find;
    
    requires 'name';
    requires '_num_steps';
    requires 'description';
    requires 'steps';
    
    our %pipeline_modules = map { $_ => 1 } findallmod(VRPipe::Pipelines);
    
    method num_steps {
        return $self->_num_steps;
    }
    method _increment_steps (PositiveInt $step_num) {
        return $self->_num_steps($step_num);
    }
    
    method add_step (VRPipe::Step $step) {
        my $step_num = $self->num_steps() + 1;
        my $sm = VRPipe::StepMember->get(step => $step, pipeline => $self, step_number => $step_num);
        $self->_increment_steps($step_num);
        $self->update;
        return $sm;
    }
    
    method _construct_pipeline (ArrayRef[VRPipe::Step] $steps, ArrayRef[VRPipe::StepAdaptorDefiner] $adaptor_defs, ArrayRef[VRPipe::StepBehaviourDefiner] $behaviour_defs) {
        # first check the pipeline hasn't already been constructed correctly
        my $all_ok = 1;
        my $schema =  $self->result_source->schema;
        my $step_num = 0;
        my @sms;
        foreach my $step (@$steps) {
            my $step_id = $step->id;
            my $rs = $schema->resultset('StepMember')->search(
                {
                    'step' => $step_id,
                    'pipeline' => $self->id,
                    'step_number' => ++$step_num
                }
            );
            
            my @results;
            while (my $sm = $rs->next) {
                push(@results, $sm);
            }
            
            if (@results == 1) {
                push(@sms, @results);
            }
            else {
                $all_ok = 0;
            }
        }
        
        unless ($all_ok) {
            # delete all stepmembers currently associated with the pipeline, as
            # long as all stepstates associated with those stepmembers are
            # complete (else throw)
            # *** if you alter a pipeline mid-run, are you forced to restart from scratch?!
            #*** not yet implemented
            
            # construct pipeline from scratch
            $self->_num_steps(0);
            $self->update;
            @sms = ();
            foreach my $step (@$steps) {
                push(@sms, $self->add_step($step));
            }
        }
        
        #*** how do we delete adaptors and behaviour we no longer want?
        
        # create adaptors
        foreach my $definer (@{$adaptor_defs}) {
            $definer->define($self);
        }
        
        # create behaviours
        foreach my $definer (@{$behaviour_defs}) {
            $definer->define($self);
        }
        
        return @sms;
    }
    
    before steps {
        if ($self->isa('VRPipe::Persistent')) {
            my $name = $self->name;
            my $module = "VRPipe::Pipelines::$name";
            
            if (exists $pipeline_modules{$module}) {
                eval "require $module;";
                unless ($@) {
                      my $obj = $module->new();
                      $self->_construct_pipeline($obj->_step_list);
                }
            } else {
                return;
            }
        }
    }
}

1;