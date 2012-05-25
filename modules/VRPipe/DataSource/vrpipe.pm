use VRPipe::Base;

class VRPipe::DataSource::vrpipe with VRPipe::DataSourceRole {
    use Digest::MD5 qw(md5_hex);
    
    method description {
        return "Use files created by VRPipe pipelines as a datasource.";
    }
    method source_description {
        return "List of pipelinesetup names or ids separated by a pipe character '|'. Each pipeline setup or id may be followed in square brackets by the step number or step name which produced the output files to be used. Optionally, if a step produced multiple types of output files, the step name or number may be followed by the output file key to identify the correct files to be used. Multiple steps from the same pipeline setup may be specified. e.g. pipeline_setup_id[step_name:output_file_key]|pipeline_setup_name[step_number1,step_number2]'";
    }
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will consist of the output files from the vrpipe datasource. If the maintain_element_grouping option is set to 1 (default), then all files produced by a dataelement in the source will be grouped into a dataelement. Otherwise, each source file will be it's own dataelement. The filter option is a string of the form 'metadata_key#regex'. If the filter_after_grouping option is set (the default), grouping based on metadata will be performed first and then the filter applied with it only being necessary for one file in the group to pass the filter by having metadata matching the regex. If the filter_after_grouping option is not set, only files which match the regex will be included and grouped based on their metadata.";
        }
        elsif ($method eq 'group_by_metadata') {
            return "Files from the source will be grouped according to their metadata keys. Requires the metadata_keys option which is a '|' separated list of metadata keys by which dataelements will be grouped. e.g. metadata_keys => 'sample|platform|library' will groups all elements with the same sample, platform and library into one dataelement. The filter option is a string of the form 'metadata_key#regex'. If the filter_after_grouping option is set (the default), grouping based on metadata will be performed first and then the filter applied with it only being necessary for one file in the group to pass the filter by having metadata matching the regex. If the filter_after_grouping option is not set, only files which match the regex will be included and grouped based on their metadata.";
        }
        
        return '';
    }
    
    has 'vrpipe_sources' => (is => 'ro',
                             isa => 'HashRef',
                             builder => '_build_vrpipe_sources',
                             lazy => 1);
    
    method _build_vrpipe_sources {
        my $schema = $self->_handle;
        my @sources = split /\|/, $self->source;
        
        my %vrpipe_sources;
        foreach my $source (@sources) {
            my ($pipeline_setup, $steps) = $source =~ m/^(.+)\[(.+)\]$/;
            $pipeline_setup ||= $source;
            my $setup;
            if ($pipeline_setup =~ /^\d+$/) {
                $setup = $schema->resultset("PipelineSetup")->find({ id => $pipeline_setup });
                unless ($setup) {
                    $self->throw("$pipeline_setup is not a valid pipeline setup id\n");
                }
            }
            else {
                $setup = $schema->resultset("PipelineSetup")->find({ name => $pipeline_setup });
                unless ($setup) {
                    $self->throw("$pipeline_setup is not a valid pipeline setup name\n");
                }
            }
            my $setup_id = $setup->id;
            
            # select desired step members
            my @step_members = $setup->pipeline->step_members;
            my %desired_steps;
            if ($steps) {
                my @steps = split /,/, $steps;
                foreach my $step_name (@steps) {
                    my ($name, $kind) = split /\:/, $step_name;
                    $kind ||= 'all';
                    if ($name =~ /^\d+$/) {
                        $desired_steps{numbers}->{$name}->{$kind} = 1;
                    }
                    else {
                        $desired_steps{names}->{$name}->{$kind} = 1;
                    }
                }
            }
            else {
                foreach my $stepm (@step_members) {
                    $desired_steps{names}->{$stepm->step->name}->{all} = 1;
                }
            }
            
            
            my $max_step = 0;
            my $final_smid;
            foreach my $stepm (@step_members) {
                my $smid = $stepm->id;
                my $step_name = $stepm->step->name;
                my $step_num = $stepm->step_number;
                if (exists $desired_steps{names}->{$step_name}) {
                    foreach my $kind (keys %{$desired_steps{names}->{$step_name}}) {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                    }
                    if ($max_step < $step_num) {
                        $max_step = $step_num;
                        $final_smid = $smid;
                    }
                }
                if (exists $desired_steps{numbers}->{$step_num}) {
                    foreach my $kind (keys %{$desired_steps{numbers}->{$step_num}}) {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                    }
                    if ($max_step < $step_num) {
                        $max_step = $step_num;
                        $final_smid = $smid;
                    }
                }
            }
            $vrpipe_sources{$setup_id}->{total_steps} = scalar @step_members;
            $vrpipe_sources{$setup_id}->{num_steps} = $max_step;
            $vrpipe_sources{$setup_id}->{final_smid} = $final_smid;
        }
        return \%vrpipe_sources;
    }
    
    method _open_source {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }
    
    method all (Defined :$handle!, Bool :$maintain_element_grouping = 1, Str :$filter?, Bool :$filter_after_grouping = 1) {
        my %args = (handle => $handle, maintain_element_grouping => $maintain_element_grouping, filter_after_grouping => $filter_after_grouping);
        if ($filter) {
            $args{filter} = $filter;
        }
        my @elements;
        foreach my $result (@{$self->_all_results(%args, complete_elements => 1)}) {
            if ($filter) {
                next unless $result->{pass_filter};
            }
            my $res = { paths => $result->{paths} };
            if ($maintain_element_grouping) {
                $res->{lane} = $result->{result}->{lane} if (exists $result->{result}->{lane});
                $res->{group} = $result->{result}->{group} if (exists $result->{result}->{group});
            }
            push @elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => $res, withdrawn => 0);
            VRPipe::DataElementLink->get(pipelinesetup => $result->{parent}->{setup_id}, parent => $result->{parent}->{element_id}, child => $elements[-1]->id);
        }
        return \@elements;
    }
    
    method _all_results (Defined :$handle!, Bool :$maintain_element_grouping = 1, Str :$filter?, Bool :$complete_elements = 1, Bool :$complete_all = 0, Bool :$filter_after_grouping = 1) {
        my ($key, $regex);
        if ($filter) {
            ($key, $regex) = split('#', $filter);
            $self->throw("Option 'filter' for vrpipe datasource was not properly formed\n") unless ($key && $regex);
        }
        
        my @output_files;
        my $vrpipe_sources = $self->vrpipe_sources;
        foreach my $setup_id (keys %{$vrpipe_sources}) {
            my $stepmembers = $vrpipe_sources->{$setup_id}->{stepmembers};
            my $setup = VRPipe::PipelineSetup->get(id => $setup_id);
            my $elements = $setup->datasource->elements;
            my $total_steps = $vrpipe_sources->{$setup_id}->{total_steps};
            
            my @per_element_output_files;
            ELEMENT: foreach my $element (@$elements) {
                my $element_state = VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $element);
                
                # complete_elements means we do not consider a dataelement until
                # it has completed its pipeline, even if the output file we want
                # comes from an earlier step. This is because our pipeline might
                # delete this file, whilst the parent pipeline is still running
                # and needs it to complete
                unless ($element_state->completed_steps == $total_steps) {
                    return([]) if $complete_all;
                    next if $complete_elements;
                }
                
                my $element_id = $element->id;
                my %element_hash = ( parent => { element_id => $element_id, setup_id => $setup_id }, result => $element->result );
                foreach my $smid (keys %{$stepmembers}) {
                    my $stepm = VRPipe::StepMember->get(id => $smid);
                    my $force = exists $stepmembers->{$smid}->{all};
                    
                    my $stepstate = $handle->resultset("StepState")->find({stepmember => $smid, dataelement => $element_id, pipelinesetup => $setup_id});
                    unless ($stepstate && $stepstate->complete) {
                        # this shouldn't happen since we did the completed_steps
                        # check above, but just incase we have an
                        # inconsistency...
                        return([]) if $complete_all;
                        next ELEMENT if $complete_elements;
                    }
                    
                    my $step_outs = $stepstate->output_files;
                    while (my ($kind, $files) = each %{$step_outs}) {
                        unless ($force) {
                            next unless exists $stepmembers->{$smid}->{$kind};
                        }
                        
                        foreach my $file (@$files) {
                            my $pass_filter = 0;
                            my $meta = $file->metadata;
                            if ($filter) {
                                # if "filter_after_grouping => 0", we filter before grouping
                                # by skipping files which don't match the regex or don't
                                # have the required metadata
                                if (defined $meta->{$key}) {
                                    $pass_filter = $meta->{$key} =~ m/$regex/;
                                    next if (!$filter_after_grouping && !$pass_filter);
                                } else {
                                    next unless $filter_after_grouping;
                                }
                            }
                            
                            if ($maintain_element_grouping) {
                                push @{$element_hash{paths}}, $file->path->stringify;
                                $element_hash{pass_filter} ||= $pass_filter;
                            }
                            else {
                                my %hash;
                                $hash{paths} = [ $file->path->stringify ];
                                $hash{metadata} = $meta if (keys %{$meta});
                                $hash{parent} = { element_id => $element_id, setup_id => $setup_id };
                                $hash{pass_filter} ||= $pass_filter;
                                push(@per_element_output_files, \%hash);
                            }
                        }
                    }
                }
                push(@per_element_output_files, \%element_hash) if ($maintain_element_grouping && exists $element_hash{paths});
            }
            push(@output_files, @per_element_output_files);
        }
        return \@output_files;
    }
    
    method group_by_metadata (Defined :$handle!, Str :$metadata_keys!, Str :$filter?, Bool :$filter_after_grouping = 1) {
        my %args = (handle => $handle, maintain_element_grouping => 0, filter_after_grouping => $filter_after_grouping);
        if ($filter) {
            $args{filter} = $filter;
        }
        
        my $group_hash;
        my @meta_keys = split /\|/, $metadata_keys;
        foreach my $hash_ref (@{$self->_all_results(%args, complete_all => 1)}) {
            my @group_keys;
            foreach my $key (@meta_keys) {
                $self->throw("Metadata key $key not present in file ".$hash_ref->{paths}->[0]."\n") unless (exists $hash_ref->{metadata}->{$key});
                push @group_keys, $hash_ref->{metadata}->{$key};
            }
            my $group_key = join '|', @group_keys;
            push(@{$group_hash->{$group_key}->{paths}}, @{$hash_ref->{paths}});
            push(@{$group_hash->{$group_key}->{parents}}, $hash_ref->{parent});
            $group_hash->{$group_key}->{pass_filter} ||= $hash_ref->{pass_filter};
        }
        
        my @elements;
        foreach my $group (sort keys %{$group_hash}) {
            my $hash_ref = $group_hash->{$group};
            if ($filter) {
                next unless $hash_ref->{pass_filter};
            }
            push @elements, VRPipe::DataElement->get(datasource => $self->_datasource_id, result => { paths => $hash_ref->{paths}, group => $group }, withdrawn => 0);
            foreach my $parent (@{$hash_ref->{parents}}) {
                VRPipe::DataElementLink->get(pipelinesetup => $parent->{setup_id}, parent => $parent->{element_id}, child => $elements[-1]->id);
            }
        }
        
        return \@elements;
    }
    
    # The changed marker for vrpipe datasources will be a comma separated list of the number 
    # of elements in each of the pipelinesetups that have completed and the number that have 
    # been withdrawn. Any change in these numbers will lead to the dataelements being revised. 
    method _element_state_status_checksum {
        my $schema = $self->_handle;
        my $sources = $self->vrpipe_sources;
        my @complete_list;
        foreach my $setup_id (sort keys %{$sources}) {
            my $num_steps = $sources->{$setup_id}->{total_steps};
            my $num_complete = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup_id, completed_steps => {'>=', $num_steps}, 'dataelement.withdrawn' => 0}, { join => 'dataelement' });
            my $num_withdrawn = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup_id, completed_steps => {'>=', $num_steps}, 'dataelement.withdrawn' => 1}, { join => 'dataelement' });
            push @complete_list, ($num_complete, $num_withdrawn);
        }
        my $digest = md5_hex join(',', @complete_list);
        return $digest;
    }
    
    method _all_datasource_files {
        my $schema = $self->_handle;
        my $sources = $self->vrpipe_sources;
        my @all_files;
        foreach my $setup_id (sort keys %{$sources}) {
            # select s,mtime from file join stepoutputfile on stepoutputfile.file=file.id
            # join stepstate on stepstate.id=stepoutputfile.stepstate
            # join dataelement on stepstate.dataelement=dataelement.id
            # where dataelement.withdrawn=0 and
            # stepoutputfile.output_key != "temp" and
            # stepstate.stepmember=$stepm and stepstate.pipelinesetup=$setup_id and 
            # stepstate.complete=1;
            my $stepm = $sources->{$setup_id}->{final_smid};
            my $rs = $schema->resultset('StepOutputFile')->search({ 'stepstate.complete' => 1, 
                                                                    'stepstate.pipelinesetup' => $setup_id, 
                                                                    'dataelement.withdrawn' => 0, 
                                                                    'stepstate.stepmember' => $stepm, 
                                                                    'output_key' => { '!=', 'temp' } }, 
                                                                  { join => {'stepstate' => 'dataelement'} });
            while (my $output = $rs->next) {
                push @all_files, $output->file->resolve;
            }
        }
        return \@all_files;
    }
    
    # The changed marker updates when file sizes or mtimes for 
    # source files change. Since mtime is updated when metadata
    # on a file is changed, 
    method _element_state_file_checksum {
        my $files  = $self->_all_datasource_files;
        my @sizes  = map { $_->s } @$files;
        my @mtimes = map { $_->mtime } @$files;
        my $digest = md5_hex join( map { defined $_ ? $_ : 'NULL' } (@sizes, @mtimes) );
        return $digest;
    }
    
    # ************************************
    # We use _element_state_status_checksum for the time being. 
    # Need to find a way to detect when files in the vrpipe
    # datasource have changed - including when metadata changes.
    method _has_changed {
        my $old_complete = $self->_changed_marker || return 1;
        my $new_complete = $self->_element_state_status_checksum;
        # my $new_complete = $self->_element_state_file_checksum;
        return ($new_complete ne $old_complete) ? 1 : 0;
    }
    
    method _update_changed_marker {
        my $complete = $self->_element_state_status_checksum;
        # my $complete = $self->_element_state_file_checksum;
        $self->_changed_marker($complete);
    }
}

1;