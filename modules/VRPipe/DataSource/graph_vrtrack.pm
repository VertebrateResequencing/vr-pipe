
=head1 NAME

VRPipe::DataSource::graph_vrtrack - get pipeline inputs from vrtrack in graph
db

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Use this datasource in a setup that follows another setup that used an irods
datasource and an *_with_warehouse_metadata method. Ie. The graph database
needs to have cram file nodes attached to VRTrack schema nodes representing the
sequencing hierarchy, and optionally QC-related nodes that the
npg_cram_stats_parser step adds.

This particular arrangement is Sanger-specific, but if you have your own way to
populate the graph in the same way, then this datasource can still be used.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::DataSource::graph_vrtrack with VRPipe::DataSourceFilterRole {
    use Digest::MD5 qw(md5_hex);
    
    our $schema;
    our @metadata_keys_to_check_for_changes = qw(library library_id id_run sample sample_accession_number sample_id reference md5 donor_id public_name sample_public_name sample_cohort sample_control);
    our @comps = (['sample', 'sample_name'], ['is_paired_read', 'lane_is_paired_read'], ['library', 'library_name'], ['id_run', 'lane_run'], ['study', 'study_name'], ['study_title', 'study_name']);
    
    has '_cached_files' => (
        is        => 'rw',
        isa       => 'ArrayRef',
        clearer   => '_clear_cache',
        predicate => '_cached',
    );
    
    method description {
        return "Use an the VRTrack schema in the graph database to extract information from";
    }
    
    method source_description {
        return "A description of a parent node(s) of your desired Lanes, in the form 'Label#property#value1,value2'. Eg. 'Study#id#123,456,789' to work with all lanes under studies with those 3 ids.";
    }
    
    method method_description (Str $method) {
        if ($method eq 'lanelet_crams') {
            return "An element will correspond to the cram file directly related to one of the lane nodes that is a child of the parent node(s) defined in the source. The group_by_metadata option will first group lanes together if they share parent nodes that you specify here; options are: Library, Sample, Study, Taxon, Gender and Alignment (the reference), and you can specify more than 1, separated by commas (all must be shared by members of a group). The parent_filter option is a string of the form 'Label#propery#value'; multiple filters can be separated by commas. The filter will look for an exact match to a property of a node that the file's node is descended from, eg. specify Sample#qc_failed#0 to only have files related to samples that have not been qc failed. The qc_filter option lets you filter on properties of the file itself or on properties of certain qc-related nodes that are children of the file and may be created by some downstream analysis; these are specified in the form 'psuedoLabel#property#operator#value', and multiple of these can be separated by commas. An example might be 'stats#sequences#>#10000,stats#reads QC failed#<#1000,genotype#pass#=#1,verifybamid#=#1,file#manual_qc#=#1,file#vrtrack_qc_passed#=#1' to only use cram files with more than 10000 total sequences and fewer than 1000 qc failed reads (according to the Bam_Stats node), with a genotype status of 'pass' (from the Genotype node), a 'pass' from the verify bam id process (from the Verify_Bam_ID node) and with manual_qc and vrtrack_qc_passed metadata set to 1 on the node representing the cram file itself. If you have an example cram file path, you can see all available labels, properties and values you might want to filter on using 'vrpipe-fileinfo --path /irods/path_to.cram --vrtrack_metadata'";
        }
        return '';
    }
    
    method _open_source {
        return $self->source; # the source can't really be opened
    }
    
    method _has_changed {
        my $old      = $self->_changed_marker;
        my $checksum = $self->_graph_query_checksum;
        $self->_changed_marker($checksum);
        $old || return 1;
        return ($checksum ne $old) ? 1 : 0;
    }
    
    method _update_changed_marker {
        $self->_changed_marker($self->_graph_query_checksum);
    }
    
    method _graph_query_checksum {
        # get filtered files under our source
        my $files = $self->_get_files;
        
        # we don't do any actual grouping here, just form a checksum based on
        # the ordered node ids combined with their various metadata
        my $data = '';
        foreach my $file (@$files) {
            $data .= '||' . $file->uuid;
            
            if (defined $file->{group}) {
                $data .= "." . $file->{group};
            }
            
            foreach my $meta ($file->{properties}, $file->{hierarchy}, $file->{qc_meta}) {
                $meta || next;
                foreach my $key (sort keys %$meta) {
                    my $val = $meta->{$key};
                    next unless defined $val;
                    # limit to printable ascii so md5_hex subroutine will work
                    $val =~ tr/\x20-\x7f//cd;
                    $data .= "|$key:$val|";
                }
            }
            
            $data .= '||';
        }
        
        my $digest = md5_hex $data;
        return $digest;
    }
    
    method _get_files {
        return $self->_cached_files if $self->_cached;
        
        my $source            = $self->source;
        my $options           = $self->options;
        my $group_by_metadata = $options->{group_by_metadata};
        my $parent_filter     = $options->{parent_filter};
        my $qc_filter         = $options->{qc_filter};
        
        $self->debug_log("graph_vrtrack datasource will do the graph query for parent node(s) $source\n");
        
        # get the cram nodes that match source
        my $t = time();
        my $files = $self->_get_cram_nodes($source, $group_by_metadata, $parent_filter);
        # this doesn't do any grouping, but ensures the files have that kind of
        # node to group on, and adds necessary group metadata to the file nodes
        # as $node->{group}. It also adds standard hierarchy metadata as
        # $node->{hierarchy}. It doesn't return files that fail parent_filter
        my $e = time() - $t;
        $self->debug_log("graph_vrtrack _get_cram_nodes() call took $e seconds\n");
        
        # do additional filtering if requested
        if ($qc_filter) {
            $t = time();
            $files = $self->_qc_filter($files, $qc_filter);
            # files now have an additional {qc_meta} hash key containing just
            # the qc properties we've filtered on
            $e = time() - $t;
            $self->debug_log("graph_vrtrack _qc_filter() call took $e seconds for " . scalar(@$files) . " files\n");
        }
        
        # cache the files
        $files ||= [];
        $self->_cached_files($files);
        
        return $files;
    }
    
    method _get_cram_nodes (Str $source, Maybe[Str] $group_by_metadata?, Maybe[Str] $parent_filter?) {
        # Study#id#3165
        my ($label, $property, $vals) = split(/#/, $source);
        my @vals = split(/,/, $vals);
        unless ($label && $property && @vals) {
            $self->throw("Invalid source: $source");
        }
        
        my %parent_filters;
        if ($parent_filter) {
            # Sample#qc_failed#0,Sample#created_date#1426150152
            foreach my $filter (split(/,/, $parent_filter)) {
                my ($label, $property, $value) = split(/#/, $filter);
                $label = lc($label);
                push(@{ $parent_filters{$label} }, [$property, $value]);
            }
        }
        
        my %group_on = map { lc($_) => 1 } split(/,/, $group_by_metadata) if $group_by_metadata;
        
        $schema ||= VRPipe::Schema->create("VRTrack");
        
        my @lanes;
        foreach my $val (@vals) {
            my @these_lanes;
            if ($label eq 'Lane') {
                my $lane = $schema->get('Lane', { $property => $val });
                @these_lanes = ($lane) if $lane;
            }
            else {
                my $node = $schema->get($label, { $property => $val });
                $node || next;
                @these_lanes = $node->closest('VRTrack', 'Lane', direction => 'outgoing', all => 1);
            }
            
            my $passed = 0;
            my $failed = 0;
            LANE: foreach my $lane (@these_lanes) {
                if (exists $parent_filters{lane}) {
                    foreach my $filter (@{ $parent_filters{lane} }) {
                        my ($property, $value) = @$filter;
                        my $actual_val = $lane->property($property);
                        
                        if (!defined $actual_val && !$value) {
                            # allow a desired 0 to match an unspecified node
                            # property
                            $failed++;
                            next;
                        }
                        
                        if (!defined $actual_val || "$actual_val" ne "$value") {
                            $failed++;
                            next LANE;
                        }
                    }
                }
                
                my $file;
                foreach my $lane_file ($lane->related(outgoing => { namespace => 'VRPipe', label => 'FileSystemElement' })) {
                    if ($lane_file->basename =~ /\.cram$/) {
                        $file = $lane_file;
                        last;
                    }
                }
                unless ($file) {
                    $failed++;
                    next;
                }
                
                my $hierarchy = $schema->get_sequencing_hierarchy($file);
                
                my @grouping;
                foreach my $label (sort keys %$hierarchy) {
                    my $node = $hierarchy->{$label};
                    if ($label ne 'lane' && exists $parent_filters{$label}) {
                        foreach my $filter (@{ $parent_filters{$label} }) {
                            my ($property, $value) = @$filter;
                            my $actual_val = $node->property($property);
                            
                            if (!defined $actual_val && !$value) {
                                $failed++;
                                next;
                            }
                            
                            if (!defined $actual_val || "$actual_val" ne "$value") {
                                $failed++;
                                next LANE;
                            }
                        }
                    }
                    
                    while (my ($key, $val) = each %{ $node->properties }) {
                        $file->{hierarchy}->{"${label}_$key"} = $val;
                    }
                    
                    if (exists $group_on{$label}) {
                        my $unique_val = $node->unique_property();
                        push(@grouping, "$label:$unique_val");
                    }
                }
                
                if (@grouping) {
                    $file->{group} = join(';', @grouping);
                }
                
                push(@lanes, [$lane, $file]);
                $passed++;
            }
            
            my $extra = $failed ? " ($failed lanes failed the parent filter)" : '';
            $self->debug_log("graph_vrtrack got $passed lanes$extra for source node $label#$property#$val\n");
        }
        
        my @files = ();
        LANE: foreach my $ref (sort { $a->[0]->unique cmp $b->[0]->unique } @lanes) {
            push(@files, $ref->[1]);
        }
        
        return \@files;
    }
    
    method _qc_filter (ArrayRef $input_files, Str $qc_filter) {
        # stats#raw total sequences#>#10000,stats#reads QC failed#<#100,genotype#pass#=#1,verifybamid#pass#=#1,file#target#=#1
        my %filters;
        my $had_file_filter = 0;
        foreach my $filter (split(/,/, $qc_filter)) {
            my ($type, $property, $operator, $value) = split(/#/, $filter);
            $type = lc($type);
            if ($operator eq '=' || $operator eq '==') {
                $operator = 'eq';
            }
            
            if ($type eq 'file') {
                $had_file_filter = 1;
                my @ok_files;
                foreach my $file (@$input_files) {
                    my $actual_val = $file->property($property);
                    
                    if (!defined $actual_val && !$value && $operator eq 'eq') {
                        push(@ok_files, $file);
                        next;
                    }
                    
                    next if !defined $actual_val;
                    
                    if (eval("$actual_val $operator $value")) {
                        push(@ok_files, $file);
                    }
                }
                $input_files = \@ok_files;
                next;
            }
            
            push(@{ $filters{$type} }, [$property, $operator, $value]);
        }
        
        my @files;
        if (keys %filters) {
            FILE: foreach my $file (@$input_files) {
                my $qc_nodes = $schema->file_qc_nodes($file);
                my %seen_types;
                while (my ($label, $qc_node) = each %$qc_nodes) {
                    my $type;
                    if ($label eq 'bam_stats') {
                        $type = 'stats';
                    }
                    elsif ($label eq 'verify_bam_id') {
                        $type = 'verifybamid';
                    }
                    else {
                        $type = $label;
                    }
                    
                    my $filters = defined $filters{$type} ? $filters{$type} : next;
                    $seen_types{$type} = 1;
                    
                    foreach my $filter (@$filters) {
                        my ($property, $operator, $value) = @$filter;
                        my $actual_val = $qc_node->property($property);
                        
                        if (!defined $actual_val && !$value && $operator eq 'eq') {
                            # the property isn't defined on the qc node, but since
                            # user wants a false value, we consider this a match
                            $file->{qc_meta}->{"${type}_$property"} = $actual_val;
                            next;
                        }
                        
                        next FILE if !defined $actual_val;
                        
                        unless (eval("$actual_val $operator $value")) {
                            next FILE;
                        }
                        
                        $file->{qc_meta}->{"${type}_$property"} = $actual_val;
                    }
                }
                
                next unless scalar(keys %seen_types) == scalar(keys %filters);
                
                push(@files, $file);
            }
        }
        elsif ($had_file_filter) {
            return $input_files;
        }
        
        return \@files;
    }
    
    method lanelet_crams (Defined :$handle!, Str :$group_by_metadata?, Str :$parent_filter?, Str :$qc_filter?) {
        my $did = $self->_datasource_id;
        my (@element_args, $group_hash, $group_order_i);
        foreach my $result ($self->_all_files()) {
            my $protocol = $result->{protocol};
            
            if ($group_by_metadata) {
                my $group_key = $result->{group};
                push(@{ $group_hash->{$group_key}->{paths} }, @{ $result->{paths} });
                $group_hash->{$group_key}->{protocol} = $protocol if defined $protocol;
                unless (defined $group_hash->{$group_key}->{order}) {
                    $group_hash->{$group_key}->{order} = ++$group_order_i;
                }
            }
            else {
                push(@element_args, { datasource => $did, result => { paths => $result->{paths}, $protocol ? (protocol => $protocol) : () } });
            }
        }
        warn "lanelet_crams called _all_files and will now _create_elements for ", scalar(@element_args), " elements\n" if $self->debug;
        
        if ($group_by_metadata) {
            foreach my $group (sort { $group_hash->{$a}->{order} <=> $group_hash->{$b}->{order} } keys %$group_hash) {
                my $data     = $group_hash->{$group};
                my $protocol = $data->{protocol};
                push(@element_args, { datasource => $did, result => { paths => $data->{paths}, group => $group, $protocol ? (protocol => $protocol) : () } });
            }
        }
        
        $self->_create_elements(\@element_args);
    }
    
    method _all_files {
        # _get_files will get called twice in row: once to see if the datasource
        # changed, and again here; _get_files caches the result, and we clear
        # the cache after getting that data.
        my $files = $self->_get_files();
        $self->_clear_cache;
        
        my @results;
        my $anti_repeat_store = {};
        foreach my $file (@$files) {
            my $new_metadata;
            foreach my $meta ($file->{properties}, $file->{hierarchy}, $file->{qc_meta}) {
                $meta || next;
                while (my ($key, $val) = each %{$meta}) {
                    $new_metadata->{$key} = $val if defined $val;
                }
            }
            
            my $file_abs_path = delete $new_metadata->{path};
            $file_abs_path ||= $file->protocolless_path;
            my $protocol = $file->protocol;
            
            # consider type to be any if not defined in the file metadata; if
            # not a VRPipe filetype it will be treated as an any
            my $type = delete $new_metadata->{type};
            $type ||= 'any';
            
            # tweak the metadata for compatibility with certain steps
            $new_metadata->{lane}         = delete $new_metadata->{"lane_unique"};
            $new_metadata->{reads}        = delete $new_metadata->{"lane_total_reads"};
            $new_metadata->{expected_md5} = $new_metadata->{md5} if defined $new_metadata->{md5};
            
            # remove extraneous metadata
            delete $new_metadata->{uuid};
            delete $new_metadata->{basename};
            delete $new_metadata->{"gender_source_gender_md5"};
            foreach my $comp (@comps) {
                my ($keep, $delete) = @$comp;
                if (defined $new_metadata->{$keep} && defined $new_metadata->{$delete} && $new_metadata->{$keep} eq $new_metadata->{$delete}) {
                    delete $new_metadata->{$delete};
                }
            }
            
            my ($vrfile) = VRPipe::File->search({ path => $file_abs_path, $protocol ? (protocol => $protocol) : () });
            my @changed_details;
            if ($vrfile) {
                # detect metadata changes
                $vrfile = $vrfile->original;
                my $current_metadata = $vrfile->metadata;
                
                delete $new_metadata->{sample_created_date} if defined $current_metadata->{sample_created_date};
                
                my @changed;
                if ($current_metadata && keys %$current_metadata) {
                    foreach my $key (@metadata_keys_to_check_for_changes) {
                        next unless defined $current_metadata->{$key};
                        next unless defined $new_metadata->{$key};
                        if (my $diff = $self->_vals_different($current_metadata->{$key}, $new_metadata->{$key})) {
                            push(@changed_details, $diff);
                        }
                    }
                }
                
                foreach my $comp (@comps) {
                    my ($keep, $delete) = @$comp;
                    if (defined $current_metadata->{$keep} && defined $new_metadata->{$delete} && $current_metadata->{$keep} eq $new_metadata->{$delete}) {
                        delete $new_metadata->{$delete};
                    }
                }
            }
            else {
                $vrfile = VRPipe::File->create(path => $file_abs_path, type => $type, $protocol ? (protocol => $protocol) : ());
            }
            
            # if there was no metadata this will add metadata to the file.
            $vrfile->add_metadata($new_metadata, replace_data => 0);
            
            my $result_hash = { paths => [$file_abs_path], protocol => $protocol, defined $file->{group} ? (group => $file->{group}) : () };
            if (@changed_details) {
                $result_hash->{changed} = [[$vrfile, $new_metadata]];
                $self->_start_over_elements_due_to_file_metadata_change($result_hash, \@changed_details, $anti_repeat_store);
                delete $result_hash->{changed};
            }
            push(@results, $result_hash);
        }
        
        return @results;
    }
}

1;
