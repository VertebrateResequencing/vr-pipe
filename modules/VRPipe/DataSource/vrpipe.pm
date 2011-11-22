use VRPipe::Base;

class VRPipe::DataSource::vrpipe with VRPipe::DataSourceRole
{
    method description {
        return "Use files created by VRPipe pipelines as a datasource.";
    }
    method source_description {
        return "List of pipelinesetup names or ids separated by a pipe character '|'. Each pipeline setup or id may be followed in square brackets by the step number 
        or step name which produced the output files to be used. Optionally, if a step produced multiple types of output files, the step name or number may be \
        followed by the output file key to identify the correct files to be used. Multiple steps from the same pipeline setup may be specified. \
        e.g. pipeline_setup_id[step_name:output_file_key]|pipeline_setup_name[step_number1,step_number2]'";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will consist of the output files from the vrpipe datasource. If the maintain_element_grouping option is set to 1 (default), then all \
            files produced by a dataelement in the source will be grouped into a dataelement. Otherwise, each source file will be it's own dataelement.";
        }
        elsif ($method eq 'group_by_metadata') {
            return "Files from the source will be grouped according to their metadata keys. Requires the metadata_keys option which is a '|' separated list of metadata \
            keys by which dataelements will be grouped. e.g. metadata_keys => 'sample|platform|library' will groups all elements with the same sample, platform and \
            library into one dataelement.";
        }
        
        return '';
    }

    has 'vrpipe_sources' => (is => 'ro',
                               isa => 'HashRef',
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
            my $final_smid;
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
                    if ($max_step < $step_num)
                    {
                        $max_step = $step_num;
                        $final_smid = $smid;
                    }
                }
                if (exists $desired_steps{numbers}->{$step_num})
                {
                    foreach my $kind (keys %{$desired_steps{numbers}->{$step_num}})
                    {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                    }
                    if ($max_step < $step_num)
                    {
                        $max_step = $step_num;
                        $final_smid = $smid;
                    }
                }
            }
            $vrpipe_sources{$setup_id}->{num_steps} = $max_step;
            $vrpipe_sources{$setup_id}->{final_smid} = $final_smid;
        }
        return \%vrpipe_sources;
    }

    method _open_source
    {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }

    method all (Defined :$handle!, Bool :$maintain_element_grouping = 1)
    {
        my @elements;
        foreach my $result (@{$self->_all_results(handle => $handle, maintain_element_grouping => $maintain_element_grouping)})
        {
            next if $result->{incomplete};
            my $paths = $result->{paths};
            push @elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => { paths => $paths }, withdrawn => 0);
            VRPipe::DataElementLink->get(pipelinesetup => $result->{parent}->{setup_id}, parent => $result->{parent}->{element_id}, child => $elements[-1]->id);
        }
        return \@elements;
    }

    method _all_results (Defined :$handle!, Bool :$maintain_element_grouping = 1)
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
                my %element_hash = ( parent => { element_id => $element->id, setup_id => $setup_id } );
                foreach my $smid (keys %{$stepmembers})
                {
                    my $stepm = $handle->resultset("StepMember")->find({ id => $smid });
                    my $stepstate = VRPipe::StepState->get(stepmember => $stepm, dataelement => $element, pipelinesetup => $setup);
                    # my $incomplete = $element_state->completed_steps < $stepm->step_number; # is this data element incomplete?
                    my $incomplete = !($stepstate->complete); # is this data element incomplete?
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
                                $hash{parent} = { element_id => $element->id, setup_id => $setup_id };
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

    method group_by_metadata (Defined :$handle!, Str :$metadata_keys!)
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
                VRPipe::DataElementLink->get(pipelinesetup => $parent->{setup_id}, parent => $parent->{element_id}, child => $elements[-1]->id);
            }
        }
        
        return \@elements;
    }

    # The changed marker for vrpipe datasources will be a comma separated list of the number 
    # of elements in each of the pipelinesetups that have completed and the number that have 
    # been withdrawn. Any change in these numbers will lead to the dataelements being revised. 
    method _element_state_status
    {
        my $schema = $self->_handle;
        my $sources = $self->vrpipe_sources;
        my @complete_list;
        foreach my $setup_id (sort keys %{$sources})
        {
            my $num_steps = $sources->{$setup_id}->{num_steps};
            my $num_complete = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup_id, completed_steps => {'>=', $num_steps}, 'dataelement.withdrawn' => 0}, { join => 'dataelement' });
            my $num_withdrawn = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup_id, completed_steps => {'>=', $num_steps}, 'dataelement.withdrawn' => 1}, { join => 'dataelement' });
            push @complete_list, ($num_complete, $num_withdrawn);
        }
        return join ',', @complete_list;
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