use VRPipe::Base;

class VRPipe::DataElementState extends VRPipe::Persistent {
    has 'pipelinesetup' => (is => 'rw',
                            isa => Persistent,
                            coerce => 1,
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            belongs_to => 'VRPipe::PipelineSetup');
    
    has 'dataelement' => (is => 'rw',
                          isa => Persistent,
                          coerce => 1,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1,
                          belongs_to => 'VRPipe::DataElement');
    
    has 'completed_steps' => (is => 'rw',
                              isa => IntSQL[4],
                              traits => ['VRPipe::Persistent::Attributes'],
                              default => 0);
    
    __PACKAGE__->make_persistent();
    
    method start_from_scratch (ArrayRef[PositiveInt] $step_numbers?) {
        my $do_our_steps = 0;
        my %step_numbers;
        if ($step_numbers && @$step_numbers > 0) {
            %step_numbers = map { $_ => 1 } @$step_numbers;
        }
        else {
            # by default, we'll start_over all steps that produced output files
            # unique to our own stepstates
            %step_numbers = map { $_ => 1 } $self->our_step_numbers;
        }
        
        # get all the stepstates made for our dataelement and pipeline and
        # start_over() them
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('StepState')->search({ dataelement => $self->dataelement->id, pipelinesetup => $self->pipelinesetup->id });
        my @sss; # (we must go through the search results before using the objects, because in using them we disconnect from the db, breaking the search)
        while (my $ss = $rs->next) {
            next unless exists $step_numbers{$ss->stepmember->step_number};
            push(@sss, $ss);
        }
        
        foreach my $ss (@sss) {
            $ss->start_over();
        }
        
        # each start_over() call will have set completed_steps(0) on us, but
        # we'll call it explicitly incase that ever changes
        $self->completed_steps(0);
        $self->update;
    }
    
    method our_step_numbers {
        my $schema = $self->result_source->schema;
        
        # get all our stepstates
        my $rs = $schema->resultset('StepState')->search({ dataelement => $self->dataelement->id, pipelinesetup => $self->pipelinesetup->id });
        my @step_states;
        my %ss_ids;
        while (my $ss = $rs->next) {
            push(@step_states, $ss);
            $ss_ids{$ss->id} = 1;
        }
        
        my %step_nums;
        foreach my $ss (@step_states) {
            my @ofiles = $ss->_output_files;
            my $ours = 1;
            # go through all the output files of our step state
            foreach my $sof (@ofiles) {
                next if $sof->output_key eq 'temp';
                
                # check that all other step states that output this same file
                # are also our own step states
                $rs = $schema->resultset('StepOutputFile')->search({ file => $sof->file->id });
                while (my $other_sof = $rs->next) {
                    my $other_ss = $other_sof->stepstate;
                    unless (exists $ss_ids{$other_ss->id}) {
                        $ours = 0;
                        last;
                    }
                }
                
                last unless $ours;
            }
            
            $ours || next;
            $step_nums{$ss->stepmember->step_number} = 1;
        }
        
        my @step_nums = sort { $a <=> $b } keys %step_nums;
        return @step_nums;
    }
}

1;