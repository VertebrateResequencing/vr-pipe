
=head1 NAME

VRPipe::Schema::VRPipe - schemas for VRPipe's own internal model

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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
    
    my $graph      = VRPipe::Persistent::Graph->new();
    my $config     = VRPipe::Config->new();
    my $fse_labels = $graph->_labels('VRPipe', 'FileSystemElement');
    
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
        else {
            $graph_pm = $self->get('PipelineMember', { sql_id => $desired_sm_id });
        }
        
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
        my $return_leaves = defined wantarray();
        my @cypher;
        foreach my $path (@$paths) {
            my $file       = file($path);
            my @components = $file->components();
            $self->throw("$path must be absolute") unless shift(@components) eq '';
            my $basename = pop(@components);
            
            # for speed we construct a single cypher query that uses merging to
            # either get existing dirs and file, or creates them with new uuids;
            # this bypasses the normal checking we have when using add() etc.,
            # but anything else would be too slow
            
            # root node, which is / for local disc, and something like
            # ftp://user:password@ftpserver:port/ for other protocols
            my $uuid          = $self->create_uuid();
            my $root_basename = '/';
            if ($protocol) {
                # encrypt everything past the first colon in case it contains a
                # password
                my ($pro, $text) = $protocol =~ /^([^:]+):(.*)/;
                $text ||= '';
                $text &&= $config->crypter->encrypt_hex($text);
                $root_basename = $pro . ':' . $text . $root_basename;
            }
            my %params = (root_basename => $root_basename, root_uuid => $uuid);
            my $cypher = $only_get ? "MATCH (root:$fse_labels { basename: { param }.root_basename })" : "MERGE (root:$fse_labels { basename: { param }.root_basename }) ON CREATE SET root.uuid = { param }.root_uuid ";
            
            # sub dirs
            if ($basename) {
                my @chain;
                my $previous        = 'root';
                my $dir_num         = 0;
                my $developing_path = '';
                foreach my $dir (@components) {
                    $uuid = $self->create_uuid();
                    $dir_num++;
                    $params{ $dir_num . '_basename' } = $dir;
                    $params{ $dir_num . '_uuid' }     = $uuid;
                    $developing_path .= "/$dir";
                    $params{ $dir_num . '_path' } = $developing_path;
                    push(@chain, $only_get ? "-[:contains]->(`$dir_num`:$fse_labels { basename: { param }.`${dir_num}_basename` })" : "MERGE (`$previous`)-[:contains]->(`$dir_num`:$fse_labels { basename: { param }.`${dir_num}_basename`, path: { param }.`${dir_num}_path` }) ON CREATE SET `$dir_num`.uuid = { param }.`${dir_num}_uuid`");
                    $previous = $dir_num;
                }
                $cypher .= join($only_get ? '' : ' ', @chain);
                
                # leaf file or dir
                $uuid                  = $self->create_uuid();
                $params{leaf_basename} = $basename;
                $params{leaf_uuid}     = $uuid;
                $params{leaf_path}     = "$path";
                $cypher .= ($only_get ? "-[:contains]->(leaf:$fse_labels { basename: { param }.leaf_basename })" : " MERGE (`$previous`)-[:contains]->(leaf:$fse_labels { basename: { param }.leaf_basename, path: { param }.leaf_path }) ON CREATE SET leaf.uuid = { param }.leaf_uuid") . ($return_leaves ? ' RETURN leaf' : '');
            }
            else {
                # we've only been asked for the root node
                $cypher .= ' RETURN root' if $return_leaves;
            }
            
            push(@cypher, [$cypher, { param => \%params }]);
        }
        
        if ($return_cypher) {
            return \@cypher;
        }
        
        if ($return_leaves) {
            my $data = $graph->_run_cypher(\@cypher);
            if ($data && $data->{nodes}) {
                foreach my $node (@{ $data->{nodes} }) {
                    bless $node, 'VRPipe::Schema::VRPipe::FileSystemElement';
                }
                
                if (wantarray()) {
                    return @{ $data->{nodes} };
                }
                else {
                    return $data->{nodes}->[0];
                }
            }
        }
        else {
            $graph->_run_cypher(\@cypher);
        }
    }
    
    method path_to_filesystemelement (ClassName|Object $self: Str|File $path, Str :$protocol?, Bool :$only_get = 0) {
        return $self->get_or_store_filesystem_paths(["$path"], $protocol ? (protocol => $protocol) : (), only_get => $only_get);
    }
    
    method filesystemelement_to_path (ClassName|Object $self: Object $file!, Bool $no_protocol?) {
        # we store the absolute path as a property on the node, but that could
        # theoretically become de-synced with reality; this returns and sets the
        # correct current absolute path
        $file->update_from_db; # to get the latest basename
        my @dirs = $file->related(incoming => { max_depth => 500, namespace => 'VRPipe', label => 'FileSystemElement', type => 'contains' });
        @dirs = reverse(@dirs);
        my $root = shift(@dirs);
        my $path = file('', (map { $_->basename } @dirs), $file->basename)->stringify;
        $file->add_properties({ path => $path });
        
        my $root_basename = $root->{properties}->{basename};
        $root_basename =~ s/\/$//;
        
        if ($no_protocol) {
            $root_basename = '';
        }
        elsif ($root_basename) {
            # decrypt any encrypted part of the protocol
            my ($pro, $text) = $root_basename =~ /^([^:]+):(.*)/;
            $text ||= '';
            $text &&= $config->crypter->decrypt_hex($text);
            $root_basename = $pro . ':' . $text;
        }
        
        return $root_basename . $path;
    }
    
    method filesystemelement_to_protocol (ClassName|Object $self: Object $file!, Bool $just_protocol_type?) {
        my @dirs = $file->related(incoming => { max_depth => 500, namespace => 'VRPipe', label => 'FileSystemElement', type => 'contains' });
        my $root = pop(@dirs);
        
        my $root_basename = $root->{properties}->{basename};
        $root_basename =~ s/\/$//;
        my ($pro, $text) = $root_basename =~ /^([^:]+):(.*)/;
        
        if ($just_protocol_type) {
            return $pro || 'file';
        }
        else {
            return 'file:/' unless $pro;
            
            # decrypt any encrypted part of the protocol
            $text ||= '';
            $text &&= $config->crypter->decrypt_hex($text);
            return $pro . ':' . $text;
        }
    }
    
    method filesystemelement_protocol_method (ClassName|Object $self: Object $file!, Str $method!) {
        my $fp = $file->{'_fp_obj'};
        unless ($fp) {
            my $protocol = $file->protocol;
            my ($pro) = $protocol =~ /^([^:]+):/;
            $fp = VRPipe::FileProtocol->create($pro, { path => $file->protocolless_path, protocol => $protocol });
            $file->{'_fp_obj'} = $fp;
        }
        
        return $fp->$method();
    }
    
    method move_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str :$protocol?) {
        unless (ref($source)) {
            $source = $self->path_to_filesystemelement($source, $protocol ? (protocol => $protocol) : (), only_get => 1);
            $source || return;
        }
        
        my $file = file($dest);
        $self->throw("$dest must be absolute") unless $file->is_absolute;
        my $dir = $self->path_to_filesystemelement($file->dir->stringify, $protocol ? (protocol => $protocol) : ());
        $dir->relate_to($source, 'contains', selfish => 1);
        $source->add_properties({ basename => $file->basename, path => $dest });
        
        return $source;
    }
    
    method _duplicate_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str $relation, Str :$protocol?) {
        unless (ref($source)) {
            $source = $self->path_to_filesystemelement($source, $protocol ? (protocol => $protocol) : (), only_get => 1);
            $source || return;
        }
        
        my $file = file($dest);
        $self->throw("$dest must be absolute") unless $file->is_absolute;
        my $dup = $self->path_to_filesystemelement($file->stringify, $protocol ? (protocol => $protocol) : ());
        $source->relate_to($dup, $relation);
        return $dup;
    }
    
    method symlink_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str :$protocol?) {
        return $self->_duplicate_filesystemelement($source, $dest, 'symlink', $protocol ? (protocol => $protocol) : ());
    }
    
    method copy_filesystemelement (ClassName|Object $self: Str|Object $source, Str $dest, Str :$protocol?) {
        return $self->_duplicate_filesystemelement($source, $dest, 'copy', $protocol ? (protocol => $protocol) : ());
    }
    
    method parent_filesystemelement (ClassName|Object $self: Str|Object $child, Str :$protocol?) {
        unless (ref($child)) {
            $child = $self->path_to_filesystemelement($child, $protocol ? (protocol => $protocol) : (), only_get => 1);
            $child || return;
        }
        
        my @related = $child->related(incoming => { max_depth => 500, leftmost => 1, namespace => 'VRPipe', label => 'FileSystemElement', type => 'symlink|copy' });
        if (@related == 1) {
            return $related[0];
        }
        return $child;
    }
    
    # can't get around to work, so the else will copy/paste from SchemaRole :(
    method add (Str $label!, HashRef|ArrayRef[HashRef] $properties!, HashRef :$incoming?, HashRef :$outgoing?) {
        if ($label eq 'File') {
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
