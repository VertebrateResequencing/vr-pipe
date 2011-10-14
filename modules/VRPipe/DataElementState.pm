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
    
    #*** stepoutputfile rows are still there after using start_from_scratch...
    
    method start_from_scratch (ArrayRef[PositiveInt] $step_numbers?) {
        # by default, we'll start_over all steps that produced output files in
        # the default output directory for our element
        my $do_our_steps = 0;
        my %step_numbers;
        if ($step_numbers && @$step_numbers > 0) {
            %step_numbers = map { $_ => 1 } @$step_numbers;
        }
        else {
            $do_our_steps = 1;
        }
        
        # get all the stepstates made for our dataelement and pipeline and
        # start_over() them
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('StepState')->search({ dataelement => $self->dataelement->id, pipelinesetup => $self->pipelinesetup->id });
        my @sss; # (we must go through the search results before using the objects, because in using them we disconnect from the db, breaking the search)
        while (my $ss = $rs->next) {
            if ($do_our_steps) {
                my $our_output_dir = $ss->stepmember->step(step_state => $ss)->output_root;
                
                my @ofiles = $ss->_output_files;
                my $ours = 1;
                foreach my $sof (@ofiles) {
                    next if $sof->output_key eq 'temp';
                    my $path = $sof->file->path;
                    
                    unless ($path =~ /^$our_output_dir/) {
                        $ours = 0;
                        last;
                    }
                }
                
                next unless $ours;
            }
            else {
                next unless exists $step_numbers{$ss->stepmember->step_number};
            }
            
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
}

1;