
=head1 NAME

VRPipe::Schema::VRPipe - schemas for VRPipe's own internal model

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, 2015 Genome Research Limited.

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

class VRPipe::Schema::VRPipe with VRPipe::SchemaRole {
    use Path::Class;
    use VRPipe::Persistent::Graph;
    use VRPipe::Config;
    use VRPipe::FileProtocol;
    use URI::Escape;
    
    my $graph      = VRPipe::Persistent::Graph->new();
    my $config     = VRPipe::Config->new();
    my $fse_labels = $graph->_labels('VRPipe', 'FileSystemElement');
    my (undef, $fse_label) = split(/:/, $fse_labels);
    
    method schema_definitions {
        return [{
                label   => 'StepState',
                unique  => [qw(uuid)],
                indexed => [qw(sql_id)]
            },
            {
                label    => 'Pipeline',
                unique   => [qw(name)],
                required => [qw(description)]
            },
            {
                label    => 'Step',
                unique   => [qw(name)],
                required => [qw(description)] # we don't yet model all the other stuff Steps have...
            },
            {
                label   => 'PipelineMember',  # (called StepMember in old model)
                unique  => [qw(uuid)],
                indexed => [qw(sql_id)]
            },
            {
                label   => 'Datasource',
                unique  => [qw(uuid)],
                indexed => [qw(sql_id type method source)]
            },
            {
                label          => 'Options',
                unique         => [qw(uuid)],
                allow_anything => 1,
                keep_history   => 1
            },
            {
                label   => 'DataElement',
                unique  => [qw(uuid)],
                indexed => [qw(sql_id)]  # we don't do withdrawn for now, since we don't have a way to keep in sync
            },
            {
                # don't add this directly; we have a pretend 'File' label that
                # you can add, supplying path (only), and you get back a
                # FileSystemElement for your file basename attached to ones for
                # the parent dirs up to root
                label          => 'FileSystemElement',     # it could be a file or a dir, and could represent something that does not exist, so we can't test to see which
                unique         => [qw(uuid)],
                indexed        => [qw(basename md5 path)], # full path is stored just to avoid a query to find the full path based on relationships when you're given a FileSystemElement in some other query; it isn't canonical but we hope to keep it up-to-date
                allow_anything => 1,                       # allow arbitrary metadata to be stored on files/dirs
                methods        => {
                    path              => sub { __PACKAGE__->filesystemelement_to_path(shift, shift) }, # don't trust the path property - calculate instead
                    protocolless_path => sub { __PACKAGE__->filesystemelement_to_path(shift, 1) },
                    cat_cmd  => sub { __PACKAGE__->filesystemelement_protocol_method(shift, 'cat_cmd') },
                    open     => sub { __PACKAGE__->filesystemelement_protocol_method(shift, 'open') },
                    openr    => sub { __PACKAGE__->filesystemelement_protocol_method(shift, 'openr') },
                    close    => sub { __PACKAGE__->filesystemelement_protocol_method(shift, 'close') },
                    move     => sub { __PACKAGE__->move_filesystemelement(shift,            shift); },
                    protocol => sub { __PACKAGE__->filesystemelement_to_protocol(shift,     shift); }
                }
            },
            {
                label   => 'KeyVal',
                unique  => [qw(keyval_md5)],
                indexed => [qw(key val)]
            },
            {
                label    => 'AutoIncrement',
                unique   => [qw(label)],
                optional => [qw(count)]
            },
            {
                label    => 'PipelineSetup',
                unique   => [qw(id)],                              # in the future from AutoIncrement, for now, from mysql
                indexed  => [qw(name output_root user unix_group)],
                optional => [qw(description)]                      # desired_farm, controlling_farm and active can't be kept in sync for now
            },
        ];
    }
    
    # this is a stop-gap while the graph db is optional and we're not modelling
    # everything in the graph directly from the start
    method ensure_state_hierarchy (Object $ss) {
        # for speed (sacrificing the SchemaRole checks etc.) we construct one
        # massive cypher query that will create any missing parts of the state
        # hierarchy, returning a graph stepstate node that can then be related
        # to some result node by the caller; we ensure the stepstate node is
        # in turn related through to pipelinesetup, dataelement, step, pipeline
        # etc.
        my $ss_id = $ss->id;
        
        my $graph_ss = $self->get('StepState', { sql_id => $ss->id });
        return $graph_ss if $graph_ss;
        
        # create the pipeline
        my $ps            = $ss->pipelinesetup;
        my $p             = $ps->pipeline;
        my $desired_sm_id = $ss->stepmember->id;
        my $graph_p       = $self->get('Pipeline', { name => $p->name });
        my $graph_pm;
        unless ($graph_p) {
            $p->block_until_locked;
            $graph_p = $self->get('Pipeline', { name => $p->name });
            if ($graph_p) {
                $graph_pm = $self->get('PipelineMember', { sql_id => $desired_sm_id });
            }
            else {
                $graph_p = $self->add('Pipeline', { name => $p->name, description => $p->description });
                
                # create the pipeline members
                my $previous_pm = $graph_p;
                foreach my $sm (sort { $a->step_number <=> $b->step_number } $p->step_members) {
                    # get/create the step
                    my $step = $sm->step;
                    my $graph_step = $self->add('Step', { name => $step->name, description => $step->description });
                    
                    my $sm_id = $sm->id;
                    my $pm = $self->add('PipelineMember', { sql_id => $sm_id }, incoming => { type => 'next_step', node => $previous_pm });
                    $pm->relate_to($graph_step, 'step', replace => 1);
                    $previous_pm = $pm;
                    
                    if ($sm_id == $desired_sm_id) {
                        $graph_pm = $pm;
                    }
                }
            }
            $p->unlock;
        }
        else {
            $graph_pm = $self->get('PipelineMember', { sql_id => $desired_sm_id });
        }
        $graph_pm || $self->throw("Failed to get a PipelineMember with an sql id of $desired_sm_id (step " . $ss->stepmember->step_number . " of pipeline " . $p->name . ')');
        
        # create the datasource
        my $ds = $ps->datasource;
        my $graph_ds = $self->get('Datasource', { sql_id => $ds->id });
        unless ($graph_ds) {
            $graph_ds = $self->add('Datasource', { sql_id => $ds->id, type => $ds->type, method => $ds->method, source => $ds->source });
            $self->add('Options', $ds->options, incoming => { type => 'options', node => $graph_ds });
        }
        
        # create the dataelement
        my $de = $ss->dataelement;
        my $graph_de = $self->get('DataElement', { sql_id => $de->id });
        unless ($graph_de) {
            $graph_de = $self->add('DataElement', { sql_id => $de->id }, incoming => { type => 'element', node => $graph_ds });
            
            my %files;
            foreach my $file (@{ $de->files || [] }) {
                push(@{ $files{ $file->protocol } }, $file->protocolless_path);
            }
            
            while (my ($protocol, $paths) = each %files) {
                foreach my $graph_fse ($self->get_or_store_filesystem_paths($paths, $protocol ne 'file:/' ? (protocol => $protocol) : ())) {
                    $graph_de->relate_to($graph_fse, 'input_file');
                }
            }
            
            # we store keyvals (metadata) on KeyVal nodes and not as properties
            # of DataElement so that we can store multiple values per key, and
            # it makes searching globally by key or value quicker, because we
            # can't index arbitrary properties of a node.
            foreach my $kv ($de->keyvallist->keyvals) {
                my ($key, $val) = @$kv;
                my $md5 = $self->md5sum($key . '.' . $val);
                $self->add('KeyVal', { keyval_md5 => $md5, key => $key, val => $val }, incoming => { type => 'keyval', node => $graph_de });
            }
        }
        
        # create pipelinesetup; in the future we will create a user-friendly
        # unique auto-incrementing id by doing something like:
        # MERGE (id:AutoIncrement { label: 'PipelineSetup' }) ON CREATE SET id.count = 1 ON MATCH SET id.count = id.count + 1 WITH id.count AS uid, id MERGE (ps:PipelinSetup { ... }) ON CREATE SET ps.id = uid ON MATCH SET id.count = id.count - 1 RETURN ps
        # for now we'll just directly use the mysql id
        my $graph_ps = $self->get('PipelineSetup', { id => $ps->id });
        unless ($graph_ps) {
            $graph_ps = $self->add('PipelineSetup', { id => $ps->id, name => $ps->name, output_root => $ps->output_root->stringify, user => $ps->user, unix_group => $ps->unix_group });
            
            $graph_ps->relate_to($graph_p,  'pipeline',   replace => 1);
            $graph_ps->relate_to($graph_ds, 'datasource', replace => 1);
            $self->add('Options', $ps->options, incoming => { type => 'options', node => $graph_ps });
        }
        
        # finally we can create the StepState
        $graph_ss = $self->add('StepState', { sql_id => $ss->id }, incoming => { type => 'stepstate', node => $graph_ps });
        $graph_ss->relate_to($graph_pm, 'pipelinestep', replace => 1);
        $graph_ss->relate_to($graph_de, 'dataelement',  replace => 1);
        return $graph_ss;
    }
    
    # files in the graph database have unique ids that are independent of their
    # initial basename and directory. This allows them to be renamed and moved
    # without affecting our knowledge of what file was produced by what
    # stepstate, and where that file is now. protocol should be a string like
    # ftp://user:password@ftpserver:port or irods: and applies to all supplied
    # paths, indicating the files are not on the local filesystem.
    method get_or_store_filesystem_paths (ClassName|Object $self: ArrayRef[Str|File] $paths!, Str :$protocol?, Bool :$return_cypher = 0, Bool :$only_get = 0) {
        my @abs_paths;
        foreach my $path (@$paths) {
            my $file       = file($path);
            my @components = $file->components();
            $self->throw("$path must be absolute") unless shift(@components) eq '';
            push(@abs_paths, $path);
        }
        
        my $root          = $self->protocol_to_root($protocol);
        my $escaped_root  = uri_escape($root);
        my $escaped_paths = uri_escape(join('///', @abs_paths));
        
        my $db = $graph->_global_label;
        my @nodes;
        my %node_path_to_index;
        foreach my $node ($graph->_call_vrpipe_neo4j_plugin_and_parse("/get_or_store_filesystem_paths/$db/$escaped_root/$escaped_paths?only_get=$only_get", namespace => 'VRPipe', label => 'FileSystemElement')) {
            bless $node, "VRPipe::Schema::VRPipe::FileSystemElement";
            $node->{root_basename} = $protocol || 'file:/';
            my $pro = 'file';
            if ($protocol) {
                ($pro) = $protocol =~ /^([^:]+):.*/;
            }
            $node->{pro} = $pro;
            push(@nodes, $node);
            $node_path_to_index{ $node->{properties}->{path} } = $#nodes;
        }
        
        # nodes are returned via hash from plugin, so order has been lost, but
        # we should return them in the same order as the input $paths
        my @sorted_nodes;
        if (@nodes) {
            foreach my $path (@abs_paths) {
                push(@sorted_nodes, $nodes[$node_path_to_index{$path}]);
            }
        }
        
        if (wantarray()) {
            return @sorted_nodes;
        }
        else {
            return $sorted_nodes[0];
        }
    }
    
    method protocol_to_root (ClassName|Object $self: Maybe[Str] $protocol?) {
        my $encryped_protocol;
        if ($protocol && $protocol ne 'file:/') {
            # encrypt everything past the first colon in case it contains a
            # password
            my ($pro, $text) = $protocol =~ /^([^:]+):(.*)/;
            $text ||= '';
            $text &&= $config->crypter->encrypt_hex($text);
            $encryped_protocol = $pro . ':' . $text;
        }
        else {
            $encryped_protocol = 'file:';
        }
        
        # root node, which is / for local disc, and something like
        # ftp://user:password@ftpserver:port/ for other protocols
        my $root_basename = '/';
        if ($encryped_protocol ne 'file:') {
            $root_basename = $encryped_protocol . $root_basename;
        }
        
        return $root_basename;
    }
    
    method path_to_filesystemelement (ClassName|Object $self: Str|File $path, Str :$protocol?, Bool :$only_get = 0) {
        return $self->get_or_store_filesystem_paths(["$path"], $protocol ? (protocol => $protocol) : (), only_get => $only_get);
    }
    
    method filesystemelement_to_path (ClassName|Object $self: Object $file!, Bool $no_protocol?) {
        # we store the absolute path as a property on the node, but that could
        # theoretically become de-synced with reality; this returns and sets the
        # correct current absolute path
        my $db   = $graph->_global_label;
        my $fid  = $file->node_id;
        my $data = $graph->_call_vrpipe_neo4j_plugin("/filesystemelement_to_path/$fid");
        
        my $path = $data->{path};
        $file->{properties}->{path} = $path;
        
        my $root_basename = $data->{root};
        $root_basename =~ s/\/$//;
        
        if ($root_basename) {
            # decrypt any encrypted part of the protocol
            my ($pro, $text) = $root_basename =~ /^([^:]+):(.*)/;
            $text ||= '';
            $text &&= $config->crypter->decrypt_hex($text);
            $root_basename = $pro . ':' . $text;
            
            $file->{root_basename} = $root_basename;
            $file->{pro}           = $pro;
        }
        
        if ($no_protocol) {
            $root_basename = '';
        }
        
        return $root_basename . $path;
    }
    
    method filesystemelement_to_protocol (ClassName|Object $self: Object $file!, Bool $just_protocol_type?) {
        if (defined $file->{pro} && defined $file->{root_basename}) {
            return $just_protocol_type ? $file->{pro} : $file->{root_basename};
        }
        
        $self->filesystemelement_to_path($file);
        return $just_protocol_type ? $file->{pro} : $file->{root_basename};
    }
    
    method filesystemelement_protocol_method (ClassName|Object $self: Object $file!, Str $method!) {
        my $fp = $file->{'_fp_obj'};
        unless ($fp) {
            my $protocol = $file->protocol;
            my ($pro) = $protocol =~ /^([^:]+):/;
            eval {
                $fp = VRPipe::FileProtocol->create($pro, { path => $file->protocolless_path, protocol => $protocol });
                $file->{'_fp_obj'} = $fp;
            };
        }
        
        return $fp->$method();
    }
    
    sub _fse_plugin_strs {
        my ($self, $source, $protocol, $dest) = @_;
        
        my $source_str;
        my $extra        = '';
        my $root         = $self->protocol_to_root($protocol);
        my $escaped_root = uri_escape($root);
        if (ref($source)) {
            $source_str = $source->node_id();
        }
        else {
            $source_str = $source;
            $extra      = "?source_root=$escaped_root";
        }
        
        my ($dir, $basename);
        if ($dest) {
            my $destfile = file($dest);
            $self->throw("$dest must be absolute") unless $destfile->is_absolute;
            $dir      = uri_escape($destfile->dir->stringify);
            $basename = uri_escape($destfile->basename);
        }
        
        return (uri_escape($source_str), $extra, $escaped_root, $dir, $basename);
    }
    
    sub _fse_populate {
        my ($self, $node, $protocol, $source) = @_;
        
        return unless $node;
        
        bless $node, "VRPipe::Schema::VRPipe::FileSystemElement";
        $node->{root_basename} = $protocol || 'file:/';
        my $pro = 'file';
        if ($protocol) {
            ($pro) = $protocol =~ /^([^:]+):.*/;
        }
        $node->{pro} = $pro;
        
        if (ref($source)) {
            $source->{properties}->{path} = $node->{properties}->{path};
            $source->{root_basename}      = $node->{root_basename};
            $source->{pro}                = $node->{pro};
            
            if (defined $source->{'_fp_obj'}) {
                eval {
                    $node->{'_fp_obj'} = VRPipe::FileProtocol->create($pro, { path => $node->{properties}->{path}, protocol => $node->{root_basename} });
                    $source->{'_fp_obj'} = $node->{'_fp_obj'};
                };
            }
        }
    }
    
    method move_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str :$protocol?) {
        my ($escaped_source_str, $extra, $escaped_root, $dir, $basename) = $self->_fse_plugin_strs($source, $protocol, $dest);
        
        my $db = $graph->_global_label;
        my $moved = $graph->_call_vrpipe_neo4j_plugin_and_parse("/filesystemelement_move/$db/$escaped_source_str/$escaped_root/$dir/$basename$extra", namespace => 'VRPipe', label => 'FileSystemElement');
        
        $self->_fse_populate($moved, $protocol, $source);
        
        return $moved;
    }
    
    method _duplicate_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str $relation, Str :$protocol?) {
        my ($escaped_source_str, $extra, $escaped_root) = $self->_fse_plugin_strs($source, $protocol, $dest);
        $dest = uri_escape($dest);
        
        my $db = $graph->_global_label;
        my $dup = $graph->_call_vrpipe_neo4j_plugin_and_parse("/filesystemelement_duplicate/$db/$escaped_source_str/$escaped_root/$dest/$relation$extra", namespace => 'VRPipe', label => 'FileSystemElement');
        
        $self->_fse_populate($dup, $protocol);
        
        return $dup;
    }
    
    method symlink_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str :$protocol?) {
        return $self->_duplicate_filesystemelement($source, $dest, 'symlink', $protocol ? (protocol => $protocol) : ());
    }
    
    method copy_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str :$protocol?) {
        return $self->_duplicate_filesystemelement($source, $dest, 'copy', $protocol ? (protocol => $protocol) : ());
    }
    
    method parent_filesystemelement (ClassName|Object $self: Str|Object $child, Str :$protocol?) {
        my ($escaped_source_str, $extra) = $self->_fse_plugin_strs($child, $protocol);
        
        my $db = $graph->_global_label;
        my $parent = $graph->_call_vrpipe_neo4j_plugin_and_parse("/filesystemelement_parent/$db/$escaped_source_str$extra", namespace => 'VRPipe', label => 'FileSystemElement');
        return unless $parent;
        
        $self->_fse_populate($parent, $protocol);
        
        if (ref($child) && $child->{id} == $parent->{id}) {
            return $child;
        }
        return $parent;
    }
    
    # can't get around to work, so the else will copy/paste from SchemaRole :(
    method add (Str $label!, HashRef|ArrayRef[HashRef] $properties!, HashRef :$incoming?, HashRef :$outgoing?) {
        if ($label eq 'File') {
            if ($incoming || $outgoing) {
                $self->throw("incoming and outgoing not currently implemented when creating File nodes");
                #*** non-trivial to implement...
            }
            
            my (@paths, $protocol);
            foreach my $prop (ref($properties) eq 'ARRAY' ? @$properties : ($properties)) {
                push(@paths, $prop->{path});
                $protocol = $prop->{protocol} if defined $prop->{protocol};
            }
            
            return $self->get_or_store_filesystem_paths(\@paths, $protocol ? (protocol => $protocol) : ());
        }
        else {
            my $history_props;
            my @nodes = $self->_get_and_bless_nodes($label, 'add_nodes', $properties, { $incoming ? (incoming => $incoming) : (), $outgoing ? (outgoing => $outgoing) : () });
            
            if ($self->_is_historical($label)) {
                foreach my $node (@nodes) {
                    $node->_maintain_property_history(0);
                }
            }
            
            if (wantarray()) {
                return @nodes;
            }
            else {
                return $nodes[0];
            }
        }
    }
    
    # can't get around to work, so the else will copy/paste from SchemaRole :(
    method get (Str $label!, HashRef $properties?) {
        if ($label eq 'File') {
            my (@paths, $protocol);
            foreach my $prop (ref($properties) eq 'ARRAY' ? @$properties : ($properties)) {
                push(@paths, $prop->{path});
                $protocol = $prop->{protocol} if defined $prop->{protocol};
            }
            return $self->get_or_store_filesystem_paths(\@paths, $protocol ? (protocol => $protocol) : (), only_get => 1);
        }
        else {
            return $self->_get_and_bless_nodes($label, 'get_nodes', $properties ? ($properties) : ());
        }
    }
}

1;
