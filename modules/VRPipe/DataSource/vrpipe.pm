
=head1 NAME

VRPipe::DataSource::vrpipe - get pipeline inputs from another VRPipe pipeline

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::DataSource::vrpipe with VRPipe::DataSourceVRPipeRole {
    use Digest::MD5 qw(md5_hex);
    
    method description {
        return "Use files created by VRPipe pipelines as a datasource.";
    }
    
    method source_description {
        return "List of pipelinesetup names or ids separated by a pipe character '|'. Each pipeline setup or id may be followed in square brackets by the step number or step name which produced the output files to be used. Optionally, if a step produced multiple types of output files, the step name or number may be followed by the output file key to identify the correct files to be used. Multiple steps from the same pipeline setup may be specified. e.g. pipeline_setup_id[step_name:output_file_key]|pipeline_setup_name[step_number1,step_number2]'";
    }
    
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "Each element will consist of the output files from the vrpipe datasource. If the maintain_element_grouping option is set to 1 (default), then all files produced by a dataelement in the source will be grouped into a dataelement. Otherwise, each source file will be it's own dataelement. The filter option is a string of the form 'metadata_key#regex'; multiple filters can be separated by commas (and neither the keys nor regexs can include hashes or commas). If the filter_after_grouping option is set (the default) and maintain_element_grouping is in use, element grouping will be performed first and then the filter applied with it only being necessary for one file in the group to pass the filter. If the filter_after_grouping option is not set, only files which match the regex will be included and grouped with others from the same element. The graph_filter option is like the filter option, except instead of being a regex against metadata stored directly on the file, you look for an exact match to a property of a node that the file's node is descended from, specifying 'namespace#label#propery#value', eg. VRTrack#Sample#qc_failed#0 to only have files related to samples that have not been qc failed. Multiple of these can be specified by separating with commas, in which case all must be true, except if the same namespace,label,property is specified more than once, in which case at least one of those values must be present. So VRTrack#Sample#qc_failed#1,VRTrack#Lane#qc_status#passed,VRTrack#Lane#qc_status#investigate means only files related to samples that were qc failed and which are also for lanes that were either set to qc passed or investigate. The include_in_all_elements option takes values in the same format as source and results in each resulting element also having the output files of the setup(s) defined here.";
        }
        elsif ($method eq 'group_by_metadata') {
            return "Files from the source will be grouped according to their metadata keys. Requires the metadata_keys option which is a '|' separated list of metadata keys by which dataelements will be grouped. e.g. metadata_keys => 'sample|platform|library' will groups all elements with the same sample, platform and library into one dataelement. The filter option is a string of the form 'metadata_key#regex'; multiple filters can be separated by commas (and neither the keys nor regexs can include hashes or commas). If the filter_after_grouping option is set (the default), grouping based on metadata will be performed first and then the filter applied with it only being necessary for one file in the group to pass the filter by having metadata matching the regex. If the filter_after_grouping option is not set, only files which match the regex will be included and grouped based on their metadata. The graph_filter option is like the filter option, except instead of being a regex against metadata stored directly on the file, you look for an exact match to a property of a node that the file's node is descended from, specifying 'namespace#label#propery#value', eg. VRTrack#Sample#qc_failed#0 to only have files related to samples that have not been qc failed. Multiple of these can be specified by separating with commas, in which case all must be true, except if the same namespace,label,property is specified more than once, in which case at least one of those values must be present. So VRTrack#Sample#qc_failed#1,VRTrack#Lane#qc_status#passed,VRTrack#Lane#qc_status#investigate means only files related to samples that were qc failed and which are also for lanes that were either set to qc passed or investigate. The include_in_all_elements option takes values in the same format as source and results in each resulting element also having the output files of the setup(s) defined here.";
        }
        elsif ($method eq 'group_all') {
            return "All output files in the vrpipe datasource will be grouped into a single element. The filter option is a string of the form 'metadata_key#regex' which will select only files with metadata matching the regex; multiple filters can be separated by commas (and neither the keys nor regexs can include hashes or commas). The graph_filter option is like the filter option, except instead of being a regex against metadata stored directly on the file, you look for an exact match to a property of a node that the file's node is descended from, specifying 'namespace#label#propery#value', eg. VRTrack#Sample#qc_failed#0 to only have files related to samples that have not been qc failed. Multiple of these can be specified by separating with commas, in which case all must be true, except if the same namespace,label,property is specified more than once, in which case at least one of those values must be present. So VRTrack#Sample#qc_failed#1,VRTrack#Lane#qc_status#passed,VRTrack#Lane#qc_status#investigate means only files related to samples that were qc failed and which are also for lanes that were either set to qc passed or investigate.";
        }
        
        return '';
    }
    
    has 'vrpipe_sources' => (
        is      => 'ro',
        isa     => 'HashRef',
        builder => '_build_vrpipe_sources',
        lazy    => 1
    );
    
    has 'include_in_all_elements_sources' => (
        is      => 'ro',
        isa     => 'HashRef',
        builder => '_build_iiae_sources',
        lazy    => 1
    );
    
    method _build_vrpipe_sources (Maybe[Str] $source?) {
        my @sources = split(/\|/, $source || $self->source);
        
        my %vrpipe_sources;
        foreach my $source (@sources) {
            my ($pipeline_setup, $steps) = $source =~ m/^(.+)\[(.+)\]$/;
            $pipeline_setup ||= $source;
            my $setup;
            if ($pipeline_setup =~ /^\d+$/) {
                ($setup) = VRPipe::PipelineSetup->search({ id => $pipeline_setup }, { rows => 1 });
                unless ($setup) {
                    $self->throw("$pipeline_setup is not a valid pipeline setup id\n");
                }
            }
            else {
                ($setup) = VRPipe::PipelineSetup->search({ name => $pipeline_setup }, { rows => 1 });
                unless ($setup) {
                    $self->throw("$pipeline_setup is not a valid pipeline setup name\n");
                }
            }
            my $setup_id = $setup->id;
            
            # select desired step members
            my @step_members = $setup->pipeline->step_members;
            my %desired_steps;
            if (defined $steps) {
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
                    $desired_steps{names}->{ $stepm->step->name }->{all} = 1;
                }
            }
            
            my $max_step   = 0;
            my $final_smid = 0;
            if (exists $desired_steps{numbers}->{'0'}) {
                foreach my $kind (keys %{ $desired_steps{numbers}->{'0'} }) {
                    $vrpipe_sources{$setup_id}->{stepmembers}->{'0'}->{$kind} = 1;
                }
            }
            foreach my $stepm (@step_members) {
                my $smid      = $stepm->id;
                my $step_name = $stepm->step->name;
                my $step_num  = $stepm->step_number;
                if (exists $desired_steps{names}->{$step_name}) {
                    foreach my $kind (keys %{ $desired_steps{names}->{$step_name} }) {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                        $vrpipe_sources{$setup_id}->{stepmember_objects}->{$smid} = $stepm;
                    }
                    if ($max_step < $step_num) {
                        $max_step   = $step_num;
                        $final_smid = $smid;
                    }
                }
                if (exists $desired_steps{numbers}->{$step_num}) {
                    foreach my $kind (keys %{ $desired_steps{numbers}->{$step_num} }) {
                        $vrpipe_sources{$setup_id}->{stepmembers}->{$smid}->{$kind} = 1;
                        $vrpipe_sources{$setup_id}->{stepmember_objects}->{$smid} = $stepm;
                    }
                    if ($max_step < $step_num) {
                        $max_step   = $step_num;
                        $final_smid = $smid;
                    }
                }
            }
            $vrpipe_sources{$setup_id}->{total_steps} = scalar @step_members;
            $vrpipe_sources{$setup_id}->{num_steps}   = $max_step;
            $vrpipe_sources{$setup_id}->{final_smid}  = $final_smid;
            $vrpipe_sources{$setup_id}->{setup_obj}   = $setup;
        }
        return \%vrpipe_sources;
    }
    
    method _build_iiae_sources {
        my $iiae = $self->options->{include_in_all_elements};
        if ($iiae) {
            return $self->_build_vrpipe_sources($iiae);
        }
        else {
            return {};
        }
    }
    
    method _open_source {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }
    
    method all (Defined :$handle!, Bool :$maintain_element_grouping = 1, Str :$filter?, Str :$graph_filter?, Str :$include_in_all_elements?, Bool :$filter_after_grouping = 1) {
        my %args = (handle => $handle, maintain_element_grouping => $maintain_element_grouping, filter_after_grouping => $filter_after_grouping);
        if ($filter) {
            $args{filter} = $filter;
        }
        if ($graph_filter) {
            $args{graph_filter} = $graph_filter;
        }
        
        my ($paths_for_all, $parents_for_all);
        if ($include_in_all_elements) {
            ($paths_for_all, $parents_for_all) = $self->_handle_include_in_all_elements($handle);
            unless ($paths_for_all && @$paths_for_all) {
                return;
            }
        }
        
        # create the elements
        my @element_args;
        $self->{result_to_linkargs} = {};
        foreach my $result (@{ $self->_all_results(%args, complete_elements => 1) }) {
            if ($filter || $graph_filter) {
                next unless $result->{pass_filter};
            }
            my $res = { paths => [@{ $result->{paths} }, @{ $paths_for_all || [] }] };
            if ($maintain_element_grouping) {
                $res->{lane}  = $result->{result}->{lane}  if (exists $result->{result}->{lane});
                $res->{group} = $result->{result}->{group} if (exists $result->{result}->{group});
            }
            
            # in order to later lookup which parent element created which
            # element, we convert $res to a string representing the files and
            # keyvals in the element we want to create, and use that as the
            # linking key. However we need to cope with a child class overriding
            # _create_elements() and adding keyvals to $res, which would
            # invalidate the lookup. To get around this, we also temporarily
            # store the lookup in the element_args and adjust things if
            # necessary in our own override of _create_elements(). We do have to
            # hope that no child class deletes our lookup key from the element
            # args!
            my $res_to_str = $self->_res_to_str($res);
            push @element_args, { datasource => $self->_datasource_id, result => $res, lookup => $res_to_str };
            $self->{result_to_linkargs}->{$res_to_str} = { pipelinesetup => $result->{parent}->{setup_id}, parent => $result->{parent}->{element_id} };
        }
        $self->_create_elements(\@element_args);
        
        # create corresponding dataelementlinks
        my @link_args;
        my %result_to_eid = map { $_->[0] . '|' . $_->[1] => $_->[2] } @{ VRPipe::DataElement->get_column_values(['filelist', 'keyvallist', 'id'], { datasource => $self->_datasource_id, withdrawn => 0 }) || [] };
        while (my ($res, $linkargs) = each %{ $self->{result_to_linkargs} }) {
            my $child = $result_to_eid{$res} || $self->throw("No DataElement was created for result $res?");
            push(@link_args, { %$linkargs, child => $child });
            
            foreach my $parent (@{ $parents_for_all || [] }) {
                push(@link_args, { pipelinesetup => $parent->{setup_id}, parent => $parent->{element_id}, child => $child });
            }
        }
        delete $self->{result_to_linkargs};
        
        VRPipe::DataElementLink->bulk_create_or_update(@link_args);
    }
    
    method group_all (Defined :$handle!, Str :$filter?, Str :$graph_filter?) {
        my %args = (handle => $handle, maintain_element_grouping => 0, filter_after_grouping => 0);
        if ($filter) {
            $args{filter} = $filter;
        }
        if ($graph_filter) {
            $args{graph_filter} = $graph_filter;
        }
        my (@paths, @parents);
        foreach my $result (@{ $self->_all_results(%args, complete_all => 1) }) {
            push @paths,   @{ $result->{paths} };
            push @parents, $result->{parent};
        }
        return unless @paths;
        $self->_create_elements([{ datasource => $self->_datasource_id, result => { paths => \@paths }, withdrawn => 0 }]);
        
        # create corresponding dataelementlinks
        my @children = VRPipe::DataElement->search({ datasource => $self->_datasource_id, withdrawn => 0 });
        my @link_args;
        foreach my $parent (@parents) {
            foreach my $child (@children) {
                push(@link_args, { pipelinesetup => $parent->{setup_id}, parent => $parent->{element_id}, child => $child->id });
            }
        }
        VRPipe::DataElementLink->bulk_create_or_update(@link_args);
    }
    
    method _all_results (Defined :$handle!, Bool :$maintain_element_grouping = 1, Str :$filter?, Str :$graph_filter, HashRef :$sources?, Bool :$complete_elements = 1, Bool :$complete_all = 0, Bool :$filter_after_grouping = 1) {
        my ($krs, $gfs, $vrpipe_graph_schema, $graph) = $self->_parse_filters($filter, $graph_filter);
        
        my @output_files;
        my $vrpipe_sources = $sources || $self->vrpipe_sources;
        foreach my $setup_id (keys %{$vrpipe_sources}) {
            my $stepmembers     = $vrpipe_sources->{$setup_id}->{stepmembers};
            my $stepmember_objs = $vrpipe_sources->{$setup_id}->{stepmember_objects};
            my $setup           = $vrpipe_sources->{$setup_id}->{setup_obj};
            my $total_steps     = $vrpipe_sources->{$setup_id}->{total_steps};
            
            my $elements_pager = $setup->datasource->elements;
            next unless $elements_pager;
            my %element_state_completed_steps = map { $_->[0] => $_->[1] } @{ VRPipe::DataElementState->get_column_values(['dataelement', 'completed_steps'], { pipelinesetup => $setup_id }) || [] };
            
            my @per_element_output_files;
            while (my $elements = $elements_pager->next) {
                ELEMENT: foreach my $element (@$elements) {
                    # complete_elements means we do not consider a dataelement until
                    # it has completed its pipeline, even if the output file we want
                    # comes from an earlier step. This is because our pipeline might
                    # delete this file, whilst the parent pipeline is still running
                    # and needs it to complete
                    my $element_id = $element->id;
                    unless ($element_state_completed_steps{$element_id} == $total_steps) {
                        return ([]) if $complete_all;
                        next if $complete_elements;
                    }
                    
                    my %element_hash = (parent => { element_id => $element_id, setup_id => $setup_id }, result => $element->result);
                    foreach my $smid (keys %{$stepmembers}) {
                        my $force = exists $stepmembers->{$smid}->{all};
                        my $step_outs;
                        if ($smid == 0) {
                            foreach my $file (@{ $element->files || $self->throw("data element " . $element->id . " gave a result with no paths") }) {
                                push @{ $step_outs->{ $file->type } }, $file;
                            }
                        }
                        else {
                            my $stepm = $stepmember_objs->{$smid};
                            my ($stepstate) = VRPipe::StepState->search({ stepmember => $smid, dataelement => $element_id, pipelinesetup => $setup_id }, { rows => 1 });
                            unless ($stepstate && $stepstate->complete) {
                                # this shouldn't happen since we did the
                                # completed_steps check above, but just in case
                                # we have an inconsistency...
                                return ([])  if $complete_all;
                                next ELEMENT if $complete_elements;
                            }
                            $step_outs = $stepstate->output_files;
                            # We expect some output files, otherwise, what are we trying to source??
                            if ($complete_elements && !(keys %{$step_outs})) {
                                $self->throw("No output files found for stepstate " . $stepstate->id . ", dataelement " . $element_id . " from pipelinesetup " . $setup_id);
                            }
                        
                        }
                        while (my ($kind, $files) = each %{$step_outs}) {
                            unless ($force) {
                                next unless exists $stepmembers->{$smid}->{$kind};
                            }
                            
                            FILE: foreach my $file (@$files) {
                                my $pass_filter = $self->_file_filter($file, $filter_after_grouping, $krs, $gfs, $vrpipe_graph_schema, $graph);
                                next FILE unless defined($pass_filter);
                                
                                if ($maintain_element_grouping) {
                                    push @{ $element_hash{paths} }, $file->path->stringify;
                                    $element_hash{pass_filter} ||= $pass_filter;
                                }
                                else {
                                    my %hash;
                                    $hash{paths} = [$file->path->stringify];
                                    my $meta = $file->metadata;
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
            }
            push(@output_files, @per_element_output_files);
        }
        
        # empty the cache before the next possible check
        $self->_clear_file_filter_cache();
        
        return \@output_files;
    }
    
    method _handle_include_in_all_elements (Defined $handle!) {
        my ($paths_for_all, $parents_for_all);
        foreach my $result (@{ $self->_all_results(handle => $handle, complete_all => 1, sources => $self->include_in_all_elements_sources) }) {
            push(@$paths_for_all,   @{ $result->{paths} });
            push(@$parents_for_all, $result->{parent});
        }
        return ($paths_for_all, $parents_for_all);
    }
    
    method group_by_metadata (Defined :$handle!, Str :$metadata_keys!, Str :$filter?, Str :$graph_filter, Str :$include_in_all_elements?, Bool :$filter_after_grouping = 1) {
        my %args = (handle => $handle, maintain_element_grouping => 0, filter_after_grouping => $filter_after_grouping);
        if ($filter) {
            $args{filter} = $filter;
        }
        if ($graph_filter) {
            $args{graph_filter} = $graph_filter;
        }
        
        my $group_hash;
        my @meta_keys = split /\|/, $metadata_keys;
        foreach my $hash_ref (@{ $self->_all_results(%args, complete_all => 1) }) {
            my @group_keys;
            foreach my $key (@meta_keys) {
                unless (exists $hash_ref->{metadata}->{$key}) {
                    $self->warn("Metadata key $key not present in file " . $hash_ref->{paths}->[0] . "\n");
                    next;
                }
                push @group_keys, $hash_ref->{metadata}->{$key};
            }
            my $group_key = join '|', @group_keys;
            push(@{ $group_hash->{$group_key}->{paths} },   @{ $hash_ref->{paths} });
            push(@{ $group_hash->{$group_key}->{parents} }, $hash_ref->{parent});
            $group_hash->{$group_key}->{pass_filter} ||= $hash_ref->{pass_filter};
        }
        
        if ($include_in_all_elements) {
            my ($paths_for_all, $parents_for_all) = $self->_handle_include_in_all_elements($handle);
            
            unless ($paths_for_all && @$paths_for_all) {
                return;
            }
            
            while (my ($group_key, $hash) = each %$group_hash) {
                push(@{ $hash->{paths} },   @$paths_for_all);
                push(@{ $hash->{parents} }, @$parents_for_all);
            }
        }
        
        # build up the dataelement args in the normal way for _create_elements
        my @element_args;
        foreach my $group (sort keys %{$group_hash}) {
            my $hash_ref = $group_hash->{$group};
            if ($filter) {
                next unless $hash_ref->{pass_filter};
            }
            push @element_args, { datasource => $self->_datasource_id, result => { paths => $hash_ref->{paths}, group => $group } };
        }
        $self->_create_elements(\@element_args);
        
        # we also need to make a dataelementlink for each dataelement we just
        # made
        my %group_to_eids;
        foreach my $de (VRPipe::DataElement->search({ datasource => $self->_datasource_id, withdrawn => 0 }, { prefetch => 'keyvallist' })) {
            my $group = $de->keyvallist->as_hashref->{group};
            $group_to_eids{$group}->{ $de->id } = 1;
        }
        
        my @link_args;
        foreach my $group (sort keys %{$group_hash}) {
            my $hash_ref = $group_hash->{$group};
            if ($filter) {
                next unless $hash_ref->{pass_filter};
            }
            
            my %parents = map { $_->{element_id} => $_ } @{ $hash_ref->{parents} };
            
            my $child_ids = $group_to_eids{$group} || $self->throw("No DataElement was created for group $group?");
            foreach my $child_id (keys %$child_ids) {
                foreach my $key (sort keys %parents) {
                    my $parent = $parents{$key};
                    push(@link_args, { pipelinesetup => $parent->{setup_id}, parent => $parent->{element_id}, child => $child_id });
                }
            }
        }
        VRPipe::DataElementLink->bulk_create_or_update(@link_args);
    }
    
    # The changed marker for vrpipe datasources will be a comma separated list
    # of the number of elements in each of the pipelinesetups that have
    # completed and the number that have been withdrawn. Any change in these
    # numbers will lead to the dataelements being revised. We also include total
    # number of unwithdrawn dataelements, to notice the change when a setup gets
    # new dataelements which could result in the withdrawal of our dataelements
    # (eg we depend on setups x and y, and y also depends on x, so when x
    # completes a new de we would create a new x-only de since y was incomplete,
    # and then while y runs on the new de we could not notice any change unless
    # we include this total number).
    # Finally, to detect changes to what passes a filter, we include the number
    # of files that pass the filter(s).
    method _element_state_status_checksum {
        my $vrpipe_sources        = $self->vrpipe_sources;
        my $iiae_sources          = $self->include_in_all_elements_sources;
        my $options               = $self->options;
        my $filter_after_grouping = defined $options->{filter_after_grouping} ? $options->{filter_after_grouping} : ($self->method eq 'group_all' ? 0 : 1);
        my $filter                = $options->{filter};
        my $graph_filter          = $options->{graph_filter};
        my ($krs, $gfs, $vrpipe_graph_schema, $graph) = $self->_parse_filters($filter, $graph_filter);
        
        my @complete_list;
        my $sanity = 0;
        foreach my $sources ($vrpipe_sources, $iiae_sources) {
            foreach my $setup_id (sort keys %{$sources}) {
                my $num_steps     = $sources->{$setup_id}->{total_steps};
                my $total_active  = VRPipe::DataElementState->search({ pipelinesetup => $setup_id, 'dataelement.withdrawn' => 0 }, { join => 'dataelement' });
                my $num_complete  = VRPipe::DataElementState->search({ pipelinesetup => $setup_id, completed_steps => { '>=', $num_steps }, 'dataelement.withdrawn' => 0 }, { join => 'dataelement' });
                my $num_withdrawn = VRPipe::DataElementState->search({ pipelinesetup => $setup_id, completed_steps => { '>=', $num_steps }, 'dataelement.withdrawn' => 1 }, { join => 'dataelement' });
                $sanity += $num_complete;
                
                my $pass_filter = -1;
                if ($filter || $graph_filter) {
                    #*** need a very quick way of checking number of files that
                    # pass the filter and graph_filter; the below is probably
                    # too slow?
                    my $files = $self->_all_datasource_files([$setup_id], $sources);
                    
                    foreach my $file (@$files) {
                        my $pass = $self->_file_filter($file, $filter_after_grouping, $krs, $gfs, $vrpipe_graph_schema, $graph);
                        $pass_filter += $pass if $pass;
                    }
                }
                
                push @complete_list, ($total_active, $num_complete, $num_withdrawn, $pass_filter);
            }
        }
        
        # it's possible for things to go wonky and for our current md5 to be
        # up to date, yet we have 0 elements; avoid the md5s matching so that
        # we'll do a full update again
        my $current = VRPipe::DataElement->search({ datasource => $self->_datasource_id, withdrawn => 0 });
        if ($sanity > 0 && $current == 0) {
            push(@complete_list, 1 + int(rand(100000 - 1)));
        }
        
        my $digest = md5_hex join(',', @complete_list);
        return $digest;
    }
    
    method _all_datasource_files (Maybe[ArrayRef[Int]] $setup_ids?, Maybe[HashRef] $sources?) {
        $sources   ||= $self->vrpipe_sources;
        $setup_ids ||= [sort keys %{$sources}];
        
        my @all_files;
        foreach my $setup_id (@$setup_ids) {
            my $pager = VRPipe::StepOutputFile->search_paged({
                    'stepstate.complete'      => 1,
                    'stepstate.pipelinesetup' => $setup_id,
                    'dataelement.withdrawn'   => 0,
                    'stepstate.stepmember'    => $sources->{$setup_id}->{final_smid},
                    'output_key'              => { '!=', 'temp' }
                },
                { join => { 'stepstate' => 'dataelement' }, prefetch => 'file' }
            );
            while (my $outputs = $pager->next) {
                foreach my $output (@$outputs) {
                    push @all_files, $output->file;
                }
            }
        }
        
        return \@all_files;
    }
    
    #*** we could use something like this instead of
    # _element_state_status_checksum, but instead of mtime we'd need to set
    # some new column value every time metadata was changed... but I don't know
    # how we'd know about changes affecting graph_filter...
    method _element_state_file_checksum {
        my $files  = $self->_all_datasource_files;
        my @sizes  = map { $_->s } @$files;
        my @mtimes = map { $_->mtime } @$files;
        my $digest = md5_hex join(map { defined $_ ? $_ : 'NULL' } (@sizes, @mtimes));
        return $digest;
    }
    
    method _has_changed {
        my $old_complete = $self->_changed_marker;
        my $new_complete = $self->_element_state_status_checksum;
        $self->_changed_marker($new_complete);
        $old_complete || return 1;
        return ($new_complete ne $old_complete) ? 1 : 0;
    }
    
    method _update_changed_marker {
        my $complete = $self->_element_state_status_checksum;
        $self->_changed_marker($complete);
    }
}

1;
