
=head1 NAME

VRPipe::DataSourceFilterRole - role for datasources that have filter options

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This role has the methods for dealing with filter options.

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

role VRPipe::DataSourceFilterRole with VRPipe::DataSourceRole {
    use VRPipe::Schema;
    
    my %file_filter_cache;
    
    method _parse_filters (Maybe[Str] $filter?, Maybe[Str] $graph_filter?) {
        my @krs;
        if ($filter) {
            foreach my $kr (split(',', $filter)) {
                my ($key, $regex) = split('#', $kr);
                $self->throw("Option 'filter' for the datasource was not properly formed\n") unless ($key && $regex);
                push(@krs, [$key, $regex]);
            }
            $self->throw("Option 'filter' for the datasource was not properly formed\n") unless @krs;
        }
        my (@gfs, $vrpipe_graph_schema, $graph);
        if ($graph_filter) {
            foreach my $gf (split(',', $graph_filter)) {
                my ($namespace, $label, $prop, $value) = split('#', $gf);
                $self->throw("Option 'graph_filter' for the datasource was not properly formed\n") unless ($namespace && $label && $prop && defined($value));
                push(@gfs, [$namespace, $label, $prop, $value]);
            }
            $self->throw("Option 'graph_filter' for the datasource was not properly formed\n") unless @gfs;
            $vrpipe_graph_schema = VRPipe::Schema->create('VRPipe');
            $graph               = $vrpipe_graph_schema->graph();
        }
        return (\@krs, \@gfs, $vrpipe_graph_schema, $graph);
    }
    
    method _file_filter (VRPipe::File|Str $file, Bool $filter_after_grouping, Maybe[ArrayRef] $krs?, Maybe[ArrayRef] $gfs?, $vrpipe_graph_schema?, $graph?) {
        # to avoid doing this twice, once while calculating the changed marker
        # and again while creating elements, we cache our results and provide
        # a _clear_file_filter_cache() method that the element creating method
        # can call before it returns. Must supply a VRPipe::File object if
        # using $krs and metadata filter, otherwise just the path as a string
        # is fine
        my $file_id = ref($file) ? $file->id : $file;
        if (exists $file_filter_cache{$file_id}) {
            return $file_filter_cache{$file_id};
        }
        
        my $meta = ref($file) ? $file->metadata : {};
        my $pass_filter = 0;
        if ($krs && @$krs) {
            # if "filter_after_grouping => 0", we filter before grouping
            # by skipping files which don't match the regex or don't
            # have the required metadata
            my $passes = 0;
            foreach my $kr (@$krs) {
                my ($key, $regex) = @$kr;
                if (defined $meta->{$key}) {
                    my $this_passed = $meta->{$key} =~ m/$regex/ ? 1 : 0;
                    if (!$filter_after_grouping && !$this_passed) {
                        $file_filter_cache{$file_id} = undef;
                        return;
                    }
                    $passes += $this_passed;
                }
                else {
                    $file_filter_cache{$file_id} = undef;
                    return unless $filter_after_grouping;
                }
            }
            $pass_filter = $passes == @$krs ? 1 : 0;
        }
        
        if ($gfs && @$gfs && (($krs && @$krs) ? $pass_filter : 1)) {
            my $file_node = $vrpipe_graph_schema->get('File', { path => ref($file) ? $file->path->stringify : $file });
            unless ($file_node) {
                $file_filter_cache{$file_id} = undef;
                return unless $filter_after_grouping;
                return 0;
            }
            
            my $passes = 0;
            foreach my $gf (@$gfs) {
                my ($namespace, $label, $prop, $value) = @$gf;
                my @nodes = $graph->related_nodes(
                    $file_node,
                    incoming => {
                        namespace => $namespace,
                        label     => $label,
                        max_depth => 20,
                        $value ? (properties => { $prop => $value }) : ()
                    }
                );
                #*** there are definitely optimisations
                # that can be made here: we can get all file nodes at once, and
                # we can check all property#value pairs on the same label at
                # once.
                if (@nodes) {
                    unless ($value) {
                        # we didn't restrict nodes to those that have $prop set
                        # so that we can now test all the nodes to see if they
                        # either have $prop set to 0, or don't have $prop set at
                        # all
                        my $ok = 0;
                        foreach my $node (@nodes) {
                            $ok++ unless $graph->node_property($node, $prop);
                        }
                        
                        # we reverse the normal logic and say all nodes must
                        # have a false/unset value, instead of only 1 node
                        # having a true value
                        if ($ok == @nodes) {
                            $passes++;
                        }
                        else {
                            $file_filter_cache{$file_id} = undef;
                            return unless $filter_after_grouping;
                        }
                    }
                    else {
                        $passes++;
                    }
                }
                else {
                    $file_filter_cache{$file_id} = undef;
                    return unless $filter_after_grouping;
                }
            }
            $pass_filter = $passes == @$gfs ? 1 : 0;
        }
        
        $file_filter_cache{$file_id} = $pass_filter;
        return $pass_filter;
    }
    
    method _clear_file_filter_cache {
        %file_filter_cache = ();
    }
}

1;
