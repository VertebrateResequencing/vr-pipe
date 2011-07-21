use VRPipe::Base;

class VRPipe::StepState extends VRPipe::Persistent {
    has 'stepmember' => (is => 'rw',
                         isa => Persistent,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1,
                         belongs_to => 'VRPipe::StepMember');
    
    has 'dataelement' => (is => 'rw',
                          isa => Persistent,
                          coerce => 1,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1,
                          belongs_to => 'VRPipe::DataElement');
    
    has 'pipelinesetup' => (is => 'rw',
                            isa => Persistent,
                            coerce => 1,
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            belongs_to => 'VRPipe::PipelineSetup');
    
    has 'output_files' => (is => 'rw',
                           isa => PersistentFileHashRef,
                           traits => ['VRPipe::Persistent::Attributes'],
                           default => sub { {} });
    
    has 'temp_files' => (is => 'rw',
                         isa => ArrayRefOfPersistent,
                         traits => ['VRPipe::Persistent::Attributes'],
                         default => sub { [] });
    
    has 'cmd_summary' => (is => 'rw',
                          isa => Persistent,
                          coerce => 1,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1,
                          belongs_to => 'VRPipe::StepCmdSummary');
    
    has 'complete' => (is => 'rw',
                       isa => 'Bool',
                       traits => ['VRPipe::Persistent::Attributes'],
                       default => 0);
    
    __PACKAGE__->make_persistent(has_many => [submissions => 'VRPipe::Submission']);
    
    method update_output_file_stats {
        my $outputs = $self->output_files;
        if ($outputs) {
            foreach my $val (values %$outputs) {
                foreach my $file (@$val) {
                    $file->update_stats_from_disc(retries => 3);
                }
            }
        }
    }
    
    method unlink_temp_files {
        foreach my $vrfile (@{$self->temp_files}) {
            $vrfile->unlink;
        }
    }
}

1;