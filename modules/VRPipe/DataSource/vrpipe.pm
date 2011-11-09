use VRPipe::Base;

# my $ds = VRPipe::DataSource->get(type => 'vrpipe',
#                                  method => 'group_by_metadata',
#                                  source => '1000 genomes phase 2 (re)mapping[17:bam,18:txt]|1000 genomes phase 2 exome mapping[17]|1000 genomes phase 2 LS454 mapping[18]',
#                                  options => {  metadata_keys => 'analysis_group|sample|platform|library' });

# source => 'pipeline_setup_name[step_num:file_type,step_num:file_type]|pipeline_setup_id[step_name]|pipeline_step_id[step_name:file_type]',

class VRPipe::DataSource::vrpipe with VRPipe::DataSourceRole
{
    has 'vrpipe_sources' => (is => 'ro',
                               isa => 'HashRef',
                               # lazy => 1,
                               builder => '_build_vrpipe_sources');

    method _build_vrpipe_sources
    {
        my $schema = $self->_handle;
        my @sources = split /\|/, $self->source;
        
        my %vrpipe_sources;
        foreach my $source (@sources)
        {
            my ($pipeline_setup, $steps) = $source =~ m/^(.+)\[(.+)\]$/;
            my $setup;
            if ($pipeline_setup =~ /^\d+$/)
            {
                $setup = $schema->resultset("PipelineSetup")->find({ id => $pipeline_setup });
                unless ($setup)
                {
                    $self->throw("$pipeline_setup is not a valid pipeline setup id\n");
                }
            }
            else {
                $setup = $schema->resultset("PipelineSetup")->find({ name => $pipeline_setup });
                unless ($setup)
                {
                    $self->throw("$pipeline_setup is not a valid pipeline setup name\n");
                }
            }
            my $setup_id = $setup->id;
            
            # select desired step members
            my %desired_steps;
            my @steps = split /,/, $steps;
            foreach my $step_name (@steps)
            {
                my ($name, $kind) = split /\:/, $step_name;
                $kind ||= 'all';
                if ($name =~ /^\d+$/)
                {
                    $desired_steps{numbers}->{$name}->{$kind} = 1;
                }
                else
                {
                    $desired_steps{names}->{$name}->{$kind} = 1;
                }
            }
            
            my $max_step = 0;
            my @step_members = $setup->pipeline->step_members;
            foreach my $stepm (@step_members)
            {
                my $smid = $stepm->id;
                my $step_name = $stepm->step->name;
                my $step_num = $stepm->step_number;
                if (exists $desired_steps{names}->{$step_name})
                {
                    foreach my $kind (keys %{$desired_steps{names}->{$step_name}})
                    {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                    }
                    $max_step = ($max_step > $step_num) ? $max_step : $step_num;
                }
                if (exists $desired_steps{numbers}->{$step_num})
                {
                    foreach my $kind (keys %{$desired_steps{numbers}->{$step_num}})
                    {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                    }
                    $max_step = ($max_step > $step_num) ? $max_step : $step_num;
                }
            }
            $vrpipe_sources{$setup_id}->{num_steps} = $max_step;
        }
        return \%vrpipe_sources;
    }

    method _open_source
    {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }

    method all (Defined :$handle, Bool :$maintain_element_grouping = 1)
    {
        my @elements;
        foreach my $result (@{$self->_all_results(handle => $handle, maintain_element_grouping => $maintain_element_grouping)})
        {
            next if $result->{incomplete};
            my $paths = $result->{paths};
            push @elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => { paths => $paths }, withdrawn => 0);
            VRPipe::DataElementLink->get(parent => $result->{parent}, child => $elements[-1]->id);
        }
        return \@elements;
    }

    method _all_results (Defined :$handle, Bool :$maintain_element_grouping = 1)
    {
        my @output_files;
        my $vrpipe_sources = $self->vrpipe_sources;
        foreach my $setup_id (keys %{$vrpipe_sources})
        {
            my $stepmembers = $vrpipe_sources->{$setup_id}->{stepmembers};
            my $setup = $handle->resultset("PipelineSetup")->find({ id => $setup_id });
            my $elements = $setup->datasource->elements;
            
            foreach my $element (@{$elements})
            {
                my $element_state = VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $element);
                my %element_hash = ( parent => $element->id );
                foreach my $smid (keys %{$stepmembers})
                {
                    my $stepm = $handle->resultset("StepMember")->find({ id => $smid });
                    my $stepstate = VRPipe::StepState->get(stepmember => $stepm, dataelement => $element, pipelinesetup => $setup);
                    my $incomplete = $element_state->completed_steps < $stepm->step_number; # is this data element incomplete?
                    my $force = exists $stepmembers->{$smid}->{all};
                    my $step_outs = $stepstate->output_files;
                    while (my ($kind, $files) = each %{$step_outs})
                    {
                        unless ($force)
                        {
                            next unless exists $stepmembers->{$smid}->{$kind};
                        }
                        
                        foreach my $file (@$files)
                        {
                            if ($maintain_element_grouping)
                            {
                                push @{$element_hash{paths}}, $file->path->stringify;
                                $element_hash{incomplete} ||= $incomplete;
                            }
                            else
                            {
                                my %hash;
                                $hash{paths} = [ $file->path->stringify ];
                                my $meta = $file->metadata;
                                $hash{metadata} = $meta if (keys %{$meta});
                                $hash{incomplete} = $incomplete;
                                $hash{parent} = $element->id;
                                push(@output_files, \%hash);
                            }
                        }
                    }
                }
                push(@output_files, \%element_hash) if ($maintain_element_grouping && exists $element_hash{paths});
            }
        }
        return \@output_files;
    }

    method group_by_metadata (Defined :$handle, Str :$metadata_keys)
    {
        my $group_hash;
        my @meta_keys = split /\|/, $metadata_keys;
        foreach my $hash_ref (@{$self->_all_results(handle => $handle, maintain_element_grouping => 0)})
        {
            my @group_keys;
            foreach my $key (@meta_keys)
            {
                $self->throw("Metadata key $key not present in file ".$hash_ref->{paths}->[0]."\n") unless (exists $hash_ref->{metadata}->{$key});
                push @group_keys, $hash_ref->{metadata}->{$key};
            }
            my $group_key = join '|', @group_keys;
            push(@{$group_hash->{$group_key}->{paths}}, @{$hash_ref->{paths}});
            push(@{$group_hash->{$group_key}->{parents}}, $hash_ref->{parent});
            $group_hash->{$group_key}->{incomplete} ||= $hash_ref->{incomplete};
        }
        
        my @elements;
        foreach my $group (sort keys %{$group_hash})
        {
            my $hash_ref = $group_hash->{$group};
            next if $hash_ref->{incomplete};
            push @elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => { paths => $hash_ref->{paths}, group => $group }, withdrawn => 0);
            foreach my $parent (@{$hash_ref->{parents}})
            {
                VRPipe::DataElementLink->get(parent => $parent, child => $elements[-1]->id);
            }
        }
        
        return \@elements;
    }

    method _element_state_status
    {
        my $schema = $self->_handle;
        my $sources = $self->vrpipe_sources;
        my @complete_list;
        # my $complete = 0;
        foreach my $setup_id (sort keys %{$sources})
        {
            my $num_steps = $sources->{$setup_id}->{num_steps};
            my $num_complete = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup_id, completed_steps => {'>=', $num_steps}, 'dataelement.withdrawn' => 0}, { join => 'dataelement' });
            my $num_withdrawn = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup_id, completed_steps => {'>=', $num_steps}, 'dataelement.withdrawn' => 1}, { join => 'dataelement' });
            push @complete_list, ($num_complete, $num_withdrawn);
            # $complete += $num_complete;
        }
        return join ',', @complete_list;
        # return $complete;
    }

    method _has_changed
    {
        my $old_complete = $self->_changed_marker || return 1;
        my $new_complete = $self->_element_state_status;
        return ($new_complete ne $old_complete) ? 1 : 0;
    }

    method _update_changed_marker
    {
        my $complete = $self->_element_state_status;
        $self->_changed_marker($complete);
    }
}

1;