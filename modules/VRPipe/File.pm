
=head1 NAME

VRPipe::File - describe and work with files on a filesystem

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

B<VRPipe> must know about all the files that are input into it, that it works
on, and that it outputs. This class stores metadata about files (both normal
filesystem metadata like size and modify time, as well as arbitrary metadata),
and also has some methods to work with files, like C<open()> and C<move()>.

It is important to always use the methods here when doing anything with a file
that B<VRPipe> has dealt with, especially in the case of moving files that are
pipeline step outputs. For example, if a step output is moved without B<VRPipe>
knowing about it, and then a subsequent step (perhaps in a different pipeline)
needs to use that file, B<VRPipe> will be forced to rerun the pipeline and
steps necessary to regenerate the required file. However if C<move()> in this
class is used to do the move, B<VRPipe> will know where the required file is
now and use it from its new location.

To avoid doing stats on disc to discover if a file exists (which can be very
slow and expensive on high-performance filesystems), this metadata is stored in
the B<VRPipe> database, so again, C<unlink()> in this class should be used if
you need to delete a file, not Perl's unlink(), or Unix B<rm>.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2015 Genome Research Limited.

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

class VRPipe::File extends VRPipe::Persistent {
    use MooseX::Aliases;
    use VRPipe::FileType;
    use VRPipe::FileProtocol;
    use File::Copy;
    use Fcntl ':mode';
    use Cwd qw(abs_path);
    use Filesys::DfPortable;
    use VRPipe::Schema;
    use VRPipe::Config;
    
    our %file_type_map = (fastq => 'fq');
    our $config = VRPipe::Config->new();
    our $global_vrpipe_schema;
    our $global_vrtrack_schema;
    
    # *** a lot of stuff depends on getting/creating files based on the path
    #     alone. However, with MySQL at least, the search on path is case
    #     insensitive, so we can't store two different files that only differ
    #     by case!
    has 'path' => (
        is      => 'rw',
        isa     => AbsoluteFile,                      # we can't be nice and auto convert relative to absolute because alterations made by moose during construction do not affect what gets put in the db
        coerce  => 1,
        traits  => ['VRPipe::Persistent::Attributes'],
        is_key  => 1,
        handles => [qw(slurp stat lstat basename dir)]
    );
    
    has 'protocol' => (
        is                   => 'rw',
        isa                  => Varchar [102],
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        default              => 'file:/',
        allow_key_to_default => 1
    );
    
    has 'type' => (
        is      => 'rw',
        isa     => FileType,
        coerce  => 1,
        traits  => ['VRPipe::Persistent::Attributes'],
        builder => '_filetype_from_extension'
    );
    
    has 'e' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        builder => 'check_file_existence_on_disc'
    );
    
    has 's' => (
        is      => 'rw',
        isa     => IntSQL [16],
        traits  => ['VRPipe::Persistent::Attributes'],
        builder => 'check_file_size_on_disc'
    );
    
    has 'mtime' => (
        is          => 'rw',
        isa         => Datetime,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        builder     => 'check_mtime_on_disc',
        is_nullable => 1
    );
    
    has 'md5' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has '_lines' => (
        is          => 'rw',
        isa         => IntSQL [16],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'keyvallist' => (
        is          => 'rw',
        isa         => Persistent,
        traits      => ['VRPipe::Persistent::Attributes'],
        belongs_to  => 'VRPipe::KeyValList',
        is_nullable => 1,
        # handles     => { meta_value => 'get_value' } *** moose complains about this, don't know why
    );
    
    has 'moved_to' => (
        is          => 'rw',
        isa         => Persistent,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        belongs_to  => 'VRPipe::File',
        is_nullable => 1
    );
    
    has 'parent' => (
        is          => 'rw',
        isa         => Persistent,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        belongs_to  => 'VRPipe::File',
        is_nullable => 1
    );
    
    has _vrpipe_schema => (
        is      => 'ro',
        isa     => 'Object',
        lazy    => 1,
        builder => '_build_vrpipe_schema',
    );
    
    has _vrtrack_schema => (
        is      => 'ro',
        isa     => 'Object',
        lazy    => 1,
        builder => '_build_vrtrack_schema',
    );
    
    has _filetype => (
        is      => 'ro',
        isa     => 'Object',
        lazy    => 1,
        builder => '_build_filetype',
    );
    
    has _fileprotocol => (
        is      => 'ro',
        isa     => 'Object',
        lazy    => 1,
        builder => '_build_fileprotocol',
        handles => [qw(cat_cmd open close pro)]
    );
    
    has _opened_for_writing => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    method _build_vrpipe_schema {
        $global_vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        return $global_vrpipe_schema;
    }
    
    method _build_vrtrack_schema {
        $global_vrtrack_schema ||= VRPipe::Schema->create('VRTrack');
        return $global_vrtrack_schema;
    }
    
    method _build_filetype {
        return VRPipe::FileType->create($self->type, { file => $self->path });
    }
    
    method _build_fileprotocol {
        my $protocol = $self->protocol;
        my ($pro) = $protocol =~ /^([^:]+):/;
        return VRPipe::FileProtocol->create($pro, { path => $self->protocolless_path->stringify, protocol => $protocol });
    }
    
    sub _stat {
        my ($self, $path) = @_;
        $path ||= $self->path;  # optional so that we can call this without a db connection by supplying the path
        $self->throw("no path!") unless $path;
        
        # in the case of permission denied we don't want to just think the file
        # doesn't exist, so we'll throw instead
        my $ok = open(my $fh, '<', $path);
        if (!$ok && $! && $! =~ /Permission denied/i) {
            $self->throw("Could not stat $path: $!");
        }
        
        # NB: if $path is a symlink, stat returns results for the actual file
        # referenced by the symlink, not the symlink itself, which is what we
        # probably want in most cases anyway
        
        # caller may call multiple *_on_disc methods in quick succession, and we
        # want to avoid doing 3 stats in a row, so we check and note if we've
        # done a stat
        my @stat;
        if (defined $self->{stated_paths} && defined $self->{stated_paths}->{$path} && $self->{stated_paths}->{$path}->[0] >= (time() - 2)) {
            @stat = @{ $self->{stated_paths}->{$path}->[1] };
            # (we don't use stat(_) in case caller is switching back and forth
            # calling different _from_disc methods on different files)
        }
        else {
            @stat = stat($path);
            $self->{stated_paths}->{$path} = [time(), \@stat];
        }
        
        return @stat;
    }
    
    sub _clear_stat_cache {
        my ($self, $path) = @_;
        $path ||= $self->path;
        $path || return;
        defined $self->{stated_paths} || return;
        delete $self->{stated_paths}->{$path};
    }
    
    method check_file_existence_on_disc (File $path?) {
        my $has_stats = $self->_stat($path);
        return $has_stats ? 1 : 0;
    }
    
    method check_file_size_on_disc (File $path?) {
        my @stat = $self->_stat($path);
        return $stat[7] ? $stat[7] : 0;
    }
    
    method check_mtime_on_disc (File $path?) {
        my @stat  = $self->_stat($path);
        my $mtime = $stat[9] ? $stat[9] : 0;
        my $dt    = DateTime->from_epoch(epoch => $mtime);
        return $dt;
    }
    
    method _filetype_from_extension {
        my $path = $self->{_column_data}->{path};
        my ($type) = $path =~ /\.([^.]*)$/;
        $type ||= 'any';
        $type = lc($type);
        if ($type eq 'gz') {
            $path =~ s/\.gz$//;
            $type = VRPipe::File->create(path => $path)->type;
        }
        if (exists $file_type_map{$type}) {
            $type = $file_type_map{$type};
        }
        eval "require VRPipe::FileType::$type;";
        return $@ ? 'any' : $type;
    }
    
    __PACKAGE__->make_persistent();
    
    around protocol {
        # decrypt possibly encrypted protocol parts
        my $protocol = $self->$orig();
        if ($protocol ne 'file:/') {
            my ($pro, $text) = $protocol =~ /^([^:]+):(.*)/;
            $text ||= '';
            $text &&= $config->crypter->decrypt_hex($text);
            $protocol = $pro . ':' . $text;
        }
        return $protocol;
    }
    
    around path {
        # prefix the path with non-standard protocol
        my $path     = $self->$orig();
        my $protocol = $self->protocol;
        if ($protocol ne 'file:/') {
            #*** Path::Class passthrough methods other than basename don't work
            # properly on complex protocols, and ftp://... stringifies to
            # ftp:/... . Not sure what to do about that...
            $path = file($protocol . $path);
        }
        return $path;
    }
    
    method protocolless_path {
        return file($self->{_column_data}->{path});
    }
    
    # we used to have a column called metadata, which was replaced with
    # keyvallist; since it's more natural to create files with a metadata arg,
    # we'll have permanent backwards-compatibility. We also encrypt parts of
    # non-default protocols in case they contain passwords.
    around bulk_create_or_update (ClassName|Object $self: @args) {
        foreach my $args (@args) {
            $self->_convert_metadata_and_encrypt_protocol($args);
        }
        return $self->$orig(@args);
    }
    
    around _get (ClassName|Object $self: Bool $create, %args) {
        $self->_convert_metadata_and_encrypt_protocol(\%args);
        return $self->$orig($create, %args);
    }
    
    around search_rs (ClassName|Object $self: HashRef $search_args, Maybe[HashRef] $search_attributes?) {
        $self->_convert_metadata_and_encrypt_protocol($search_args);
        my @args = ($search_args);
        push(@args, $search_attributes) if $search_attributes;
        return $self->$orig(@args);
    }
    
    sub _convert_metadata_and_encrypt_protocol {
        my ($self, $args) = @_;
        
        if (exists $args->{metadata}) {
            my $meta = delete $args->{metadata};
            $args->{keyvallist} = VRPipe::KeyValList->get(hash => $meta)->id;
        }
        
        if (exists $args->{protocol}) {
            my $protocol = $args->{protocol};
            my ($pro, $text) = $protocol =~ /^([^:]+):(.*)/;
            $text ||= '';
            $text &&= $config->crypter->encrypt_hex($text);
            $args->{protocol} = $pro . ':' . $text;
        }
    }
    
    around open (OpenMode $mode, Str :$permissions?, Bool :$backwards?, Int :$retry = 0) {
        my $path = $self->protocolless_path;
        
        my $pro     = $self->pro;
        my $check_e = $pro eq 'file';
        if ($check_e && ($mode eq '<' && !$self->e)) {
            $self->update_stats_from_disc;
            # give it a second chance...
            if (!$self->e) {
                $self->throw("File '$path' does not exist, so cannot be opened for reading");
            }
        }
        
        my $type = VRPipe::FileType->create($self->type, { file => $path });
        my $rs = $type->record_separator;
        unless (defined $backwards) {
            $backwards = $type->read_backwards;
        }
        
        $self->_clear_stat_cache($path);
        $self->disconnect;
        
        my $fh;
        eval { $fh = $self->$orig($mode, { permissions => $permissions, backwards => $backwards, record_separator => $rs }); };
        
        if ($fh) {
            if (index($mode, '>') == 0) {
                $self->_opened_for_writing(1);
                $self->e($self->check_file_existence_on_disc) if $check_e;
                $self->update;
            }
        }
        else {
            if ($retry > 59 && $mode eq '<') {
                $self->throw("Failed to open '$path' after multiple retries: $!");
            }
            elsif ($mode eq '>') {
                $self->throw("Could not write to '$path': $!\n");
            }
            else {
                # we think the file exists, so sleep a second and try again
                $self->disconnect;
                sleep(1);
                $self->_clear_stat_cache($path);
                $self->e($self->check_file_existence_on_disc) if $check_e;
                $self->update;
                return $self->open($mode, defined $permissions ? (permissions => $permissions) : (), defined $backwards ? (backwards => $backwards) : (), retry => ++$retry);
            }
        }
        
        $self->disconnect;
        return $fh;
    }
    
    method openr {
        return $self->open('<');
    }
    
    method openw {
        return $self->open('>');
    }
    
    around close {
        $self->$orig() || return;
        if ($self->_opened_for_writing) {
            $self->_clear_stat_cache;
            $self->update_stats_from_disc;
            $self->_opened_for_writing(0);
        }
        return 1;
    }
    
    method last_line {
        my $fh = $self->open('<', backwards => 1);
        my $line = <$fh>;
        close($fh);
        return $line;
    }
    
    # speed critical, so sub instead of method
    sub metadata {
        my ($self, $meta, %opts) = @_;
        
        if ($meta) {
            $self->keyvallist(VRPipe::KeyValList->get(hash => $meta)->id);
        }
        
        if (defined wantarray) {
            my $keyvallist = $self->keyvallist || return {};
            my $meta = $keyvallist->as_hashref;
            
            if ($opts{include_vrtrack}) {
                my $vrtrack = $self->_vrtrack_schema;
                my $vrmeta = $vrtrack->vrtrack_metadata(path => $self->protocolless_path, protocol => $self->protocol);
                while (my ($key, $val) = each %$vrmeta) {
                    $meta->{$key} = $val;
                }
            }
            
            return $meta;
        }
    }
    
    method meta_value (Str $key) {
        my $kvl = $self->keyvallist || return;
        return $kvl->get_value($key);
    }
    
    # speed critical, so sub instead of method
    sub add_metadata {
        my ($self, $meta, @args) = @_;
        my %args = (replace_data => 1, @args);
        my $replace_data = $args{replace_data};
        
        $self->block_until_locked;
        $self->reselect_values_from_db; # block_until_locked won't have done this if it didn't have to block
        
        my $transaction = sub {
            my $final_meta = $self->metadata;
            
            while (my ($key, $val) = each %$meta) {
                # we don't always overwrite existing values
                if ($replace_data) {
                    next unless (defined $val && "$val" ne "");
                }
                else {
                    next if exists $final_meta->{$key};
                }
                $final_meta->{$key} = $val;
            }
            
            $self->metadata($final_meta);
            $self->update;
        };
        $self->do_transaction($transaction, "Failed to add_metadata for file " . $self->{_column_data}->{path});
        $self->unlock;
        
        my $resolve = $self->resolve;
        if ($resolve ne $self) {
            $resolve->add_metadata($meta, replace_data => $replace_data);
        }
    }
    
    sub merge_metadata {
        my ($self, $new_meta) = @_;
        
        $self->block_until_locked;
        $self->reselect_values_from_db;
        
        my $final_meta;
        my $transaction = sub {
            my $current_meta = $self->metadata;
            
            while (my ($key, $current) = each %$current_meta) {
                my $val;
                if (defined $new_meta->{$key}) {
                    my $new          = $new_meta->{$key};
                    my @new_vals     = ref($new) && ref($new) eq 'ARRAY' ? @{$new} : ($new);
                    my @current_vals = ref($current) && ref($current) eq 'ARRAY' ? @{$current} : ($current);
                    my %current_vals = map { $_ => 1 } @current_vals;
                    
                    my @merged_vals = @current_vals;
                    foreach my $new_val (@new_vals) {
                        next if exists $current_vals{$new_val};
                        push(@merged_vals, $new_val);
                    }
                    
                    $val = [sort @merged_vals];
                }
                else {
                    $val = $current;
                }
                
                $final_meta->{$key} = $val;
            }
            
            while (my ($key, $val) = each %$new_meta) {
                next if exists $current_meta->{$key};
                $final_meta->{$key} = $val;
            }
            
            $self->metadata($final_meta);
            $self->update;
        };
        $self->do_transaction($transaction, "Failed to merge_metadata for file " . $self->{_column_data}->{path});
        $self->unlock;
        
        my $resolve = $self->resolve;
        if ($resolve ne $self) {
            $resolve->add_metadata($final_meta, replace_data => 1);
        }
    }
    
    method touch {
        $self->path->touch;
        $self->update_stats_from_disc;
    }
    
    method remove {
        my $path = $self->path;
        
        my $worked = $path->remove;
        $self->update_stats_from_disc;
        if ($worked) {
            $self->_lines(undef);
            $self->parent(undef);
            $self->md5(undef);
            $self->update;
        }
        
        #*** removed the debugging below to triple the speed
        # my %stepstates = map { $_->stepstate->id => $_->stepstate } VRPipe::StepOutputFile->search({ file => $self->id });
        # while (my ($ss_id, $ss) = each %stepstates) {
        #     $ss->pipelinesetup->log_event("File->remove() called for StepOutputFile $path, " . ($worked ? 'and it worked' : 'but it failed'), dataelement => $ss->dataelement->id, stepstate => $ss->id);
        # }
        
        return $worked;
    }
    alias unlink => 'remove';
    alias rm     => 'remove';

=head2 move
 
 Title   : move (alias mv)
 Usage   : $obj->move($dest);
 Function: Move this file to another path. dest receives $obj's metadata, and
           $obj knows it was moved_to $dest. $dest only inherits $obj's parent
           if $obj had one. Pre-checks if there is sufficient disk space at the
           destination, and also fails if destination has less that 5% disk
           space remaining globally, to avoid breaking certain filesystems.
 Returns : n/a
 Args    : VRPipe::File destination file
           optionally, check_md5s => 1 to make sure the move was perfect by
           doing a copy, md5 check and then deletion of source

=cut
    
    method move (VRPipe::File $dest, Bool :$check_md5s = 0) {
        # have we already been moved there?
        if (!$self->e && $dest->e && $self->moved_to && $self->moved_to->id == $dest->id) {
            return 1;
        }
        
        my $sp         = $self->path;
        my $dp         = $dest->path;
        my %stepstates = map { $_->stepstate->id => $_->stepstate } VRPipe::StepOutputFile->search({ file => $self->id });
        while (my ($ss_id, $ss) = each %stepstates) {
            $ss->pipelinesetup->log_event("File->move() called for StepOutputFile $sp => $dp", dataelement => $ss->dataelement->id, stepstate => $ss->id);
        }
        
        my $success;
        if ($check_md5s) {
            $success = $self->copy($dest);
        }
        else {
            # is it safe to move? (this check requires enough disk space for a
            # copy, even though we may do a direct mv requiring no additional
            # disk space if both sp and dp are on the same filesystem... but
            # whatever)
            $self->check_destination_space($dp->dir);
            
            $self->disconnect;
            my $error = '[no error msg]';
            if (-l $sp) {
                # File::Copy::move copies symlinks across filesystem boundries
                # as the files they point to instead of copying them as
                # symlinks
                my $dst = readlink($sp);
                $success = symlink($dst, $dp);
                $error = $!;
                if ($success) {
                    unlink($sp);
                }
            }
            else {
                $success = File::Copy::move($sp, $dp);
                $error = $!;
            }
            
            $dest->update_stats_from_disc;
            unless ($success) {
                $self->update_stats_from_disc;
                $self->throw("move of $sp => $dp failed: $error");
            }
            $dest->add_metadata($self->metadata);
        }
        
        if ($success) {
            my $parent = $self->parent;
            if ($parent) {
                $dest->parent($parent);
                $dest->update;
            }
            # if this file was the destination of any symlinks, update these
            # symlinks to point to the new location
            foreach my $symlink (VRPipe::File->search({ parent => $self->id })) {
                $dest->update_symlink($symlink);
            }
            
            $self->moved_to($dest);
            $self->update;
            $self->remove; # to update stats and _lines and actually delete us
            
            # update the graph db if the source was in the graph
            $self->_vrpipe_schema->move_filesystemelement($self->path->stringify, $dest->path->stringify);
            
            return 1;
        }
        else {
            my $sp = $self->path;
            my $dp = $dest->path;
            $self->throw("move of $sp => $dp failed");
        }
    }
    alias mv => 'move';

=head2 symlink
 
 Title   : symlink
 Usage   : $obj->symlink($dest);
 Function: Create a soft symlink to this file at another path.
 Returns : n/a
 Args    : VRPipe::File destination file

=cut
    
    method symlink (VRPipe::File $dest) {
        my $sp      = $self->path;
        my $dp      = $dest->path;
        my $success = symlink($sp, $dp);
        
        # allow failures due to the symlink already existing
        unless ($success) {
            if (-e $dp && -l $dp && abs_path($dp) eq $sp) {
                my $parent = $dest->parent;
                if ($parent && $parent->id == $self->id) {
                    return;
                }
                $success = 1;
            }
        }
        
        unless ($success) {
            $self->throw("symlink of $sp => $dp failed: $!");
        }
        else {
            $dest->update_stats_from_disc;
            $dest->add_metadata($self->metadata);
            $dest->parent($self);
            $dest->update;
            
            # update the graph db if the source was in the graph
            $self->_vrpipe_schema->symlink_filesystemelement($self->path->stringify, $dest->path->stringify);
        }
    }
    
    method update_symlink (VRPipe::File $dest) {
        # if the symlink was removed from disk by a process outside of
        # vrpipe, then don't replace the non-existant
        # file with a new symlink.
        # don't check $dest->e because the destination of the
        # symlink has moved and does not exist
        unless (-l $dest->path) {
            $dest->update_stats_from_disc;
            $dest->parent(undef);
            $dest->update;
            return;
        }
        
        # otherwise update the symlink
        $dest->remove;
        $self->symlink($dest);
    }

=head2 resolve
 
 Title   : resolve
 Usage   : my $real_file = $obj->resolve;
 Function: If this file was created as a symlink (using VRPipe::File->symlink),
           returns the VRPipe::File corresponding to the real file. If this file
           was moved somewhere else, returns the VRPipe::File corresponding to
           the current location. If neither of these is true, returns itself.
 Returns : VRPipe::File object
 Args    : not_symlinks => boolean (do not resolve symlinks)

=cut
    
    method resolve (Bool :$not_symlinks = 0) {
        my $links_resolved;
        my $parent = $self->parent;
        if ($parent && !$not_symlinks) {
            $links_resolved = $parent->resolve;
        }
        else {
            $links_resolved = $self;
        }
        
        my $fully_resolved;
        my $moved_to = $links_resolved->moved_to;
        if ($moved_to) {
            $fully_resolved = $moved_to->resolve(not_symlinks => $not_symlinks);
        }
        else {
            $fully_resolved = $links_resolved;
        }
        
        return $fully_resolved;
    }

=head2 original
 
 Title   : original
 Usage   : my $original_file = $obj->original;
 Function: If this file was the result of moving another file, returns the
           VRPipe::File corresponding to the original location. Otherwise,
           returns itself. Does NOT resolve a symlink to the file it links to.
 Returns : VRPipe::File object in scalar context, list of File ids (back through
           the chain of moved files to the id of the original file, excluding
           self) in list context
 Args    : n/a

=cut
    
    method original {
        my @fids;
        my $current_id = $self->id;
        while (1) {
            my ($parent_id) = VRPipe::File->get_column_values('id', { moved_to => $current_id }, { rows => 1 });
            $parent_id || last;
            push(@fids, $parent_id);
            $current_id = $parent_id;
        }
        
        if (wantarray) {
            return @fids;
        }
        elsif (@fids) {
            return VRPipe::File->get(id => $fids[-1]);
        }
        return $self;
    }

=head2 output_by
 
 Title   : output_by
 Usage   : my @step_states = $obj->output_by;
 Function: If this file was a step output, returns the stepstates that created
           it.
 Returns : list of VRPipe::StepState objects; in scalar context returns 1 if
           any objects would have been returned (not the true count)
 Args    : Boolean, which if true, only returns a single stepstate (not in a
           list), and only if there is just 1 stepstate that created it
           (discounting those that ran the exact same command line).

=cut
    
    method output_by (Bool $single = 0) {
        # resolve first, then work backwards to get all file ids that
        # represented us in the past
        my $resolved = $self->resolve(not_symlinks => 1);
        my @fids = $resolved->original;
        unshift(@fids, $resolved->id) if $resolved->id != $self->id;
        unshift(@fids, $self->id);
        
        my $quick = 0;
        if (!wantarray && !$single) {
            $quick = 1;
        }
        
        # get the stepstates
        my @sss;
        foreach my $fid (@fids) {
            if ($quick) {
                my $found = VRPipe::StepOutputFile->search({ file => $fid });
                return 1 if $found;
            }
            else {
                # try limiting to stepstates that have submissions, since this
                # is much faster in the case 1 stepstate and sub created the
                # file, but we have thousands of other stepstates that have
                # no submissions and no same_submissions_as either (?!)
                my @these_sss = map { $_->stepstate } VRPipe::StepOutputFile->search({ file => $fid, $single ? ('stepstate.same_submissions_as' => undef, 'submissions.id' => { '!=' => undef }) : () }, { join => { stepstate => 'submissions' }, prefetch => 'stepstate' });
                if ($single && !@these_sss) {
                    # however it's legitimate that we can have a sof with no
                    # submissions
                    @these_sss = map { $_->stepstate } VRPipe::StepOutputFile->search({ file => $fid, $single ? ('stepstate.same_submissions_as' => undef) : () }, { prefetch => 'stepstate' });
                }
                push(@sss, @these_sss);
            }
        }
        return 0 if $quick;
        
        # remove dups
        my %sss = map { $_->id => $_ } @sss;
        @sss = sort { $a->id <=> $b->id } values %sss;
        
        return @sss unless $single;
        
        # if all of them share the same job, return the first one
        my %jobs;
        foreach my $ss (@sss) {
            foreach my $sub ($ss->submissions) {
                $jobs{ $sub->job->id }++;
            }
        }
        while (my ($jid, $count) = each %jobs) {
            if ($count == @sss) {
                return $sss[0];
            }
        }
        
        if (keys %jobs == 0) {
            return $sss[0];
        }
        
        return;
    }

=head2 input_to
 
 Title   : input_to
 Usage   : my @dataelements = $obj->input_to;
 Function: If this file was one of the files of any dataelements, returns them.
 Returns : list of VRPipe::DataElement objects; in scalar context returns 1 if
           any objects would have been returned (not the true count)
 Args    : n/a

=cut
    
    method input_to {
        # work backwards to get all file ids that represented us in the past;
        # we don't resolve because we can't be the input to something that used
        # a file path we moved to
        my @fids = $self->original;
        unshift(@fids, $self->id);
        
        my $quick = 0;
        if (!wantarray) {
            $quick = 1;
        }
        
        # get the filelists our files are part of
        my %fls;
        foreach my $fid (@fids) {
            my @fls = VRPipe::FileListMember->get_column_values('filelist', { file => $fid });
            foreach my $fl (@fls) {
                $fls{$fl} = 1;
            }
        }
        return 0 unless keys %fls;
        
        # get the dataelements that use those filelists
        my @des;
        foreach my $flid (keys %fls) {
            if ($quick) {
                my $found = VRPipe::DataElement->search({ filelist => $flid });
                return 1 if $found;
            }
            else {
                push(@des, VRPipe::DataElement->search({ filelist => $flid }));
            }
        }
        return 0 if $quick;
        
        # remove dups
        my %des = map { $_->id => $_ } @des;
        @des = sort { $a->id <=> $b->id } values %des;
        
        return @des;
    }

=head2 copy
 
 Title   : copy (alias cp)
 Usage   : $obj->copy($dest);
 Function: Copy this file to another path, checking md5s to make sure the copy
           was perfect. $dest receives $obk's metadata. Pre-checks if there is
           sufficient disk space at the destination, and also fails if
           destination has less that 5% disk space remaining globally, to avoid
           breaking certain filesystems.
 Returns : n/a
 Args    : VRPipe::File source file, VRPipe::File destination file

=cut
    
    method copy (VRPipe::File $dest) {
        my $sp = $self->path;
        my $dp = $dest->path;
        $self->throw("source ($sp) and destination ($dp) of a copy cannot be the same") if $sp eq $dp;
        $dest->update_stats_from_disc;
        my $d_existed = $dest->e;
        
        # has it already been copied successfully?
        if ($d_existed) {
            if (!$self->check_file_existence_on_disc) {
                # ... we'll just have to hope the copy was good
                $self->warn("The destination ($dp) exists, but the source ($sp) doesn't - will assume the copy worked previously");
                return 1;
            }
            elsif ($dest->md5 && $self->md5 && $self->md5 eq $dest->md5) {
                # copy was definitely already good
                return 1;
            }
            else {
                # we can't trust the copy; delete it and we'll copy it again
                $dest->remove;
            }
        }
        
        # is it safe to copy?
        $self->check_destination_space($dest->path->dir);
        
        $self->disconnect;
        my $success;
        my $error = '[no error msg]';
        if (-l $sp) {
            # File::Copy::copy copies symlinks across filesystem boundries
            # as the files they point to instead of copying them as
            # symlinks
            my $dst = readlink($sp);
            $success = symlink($dst, $dp);
            $error = $!;
        }
        else {
            # File::Copy::copy and similar do not preserve ownership, so we use
            # unix cp -p instead, which is --preserve=mode,ownership,timestamps
            $success = !system("cp -p $sp $dp");
            $error   = $!;
        }
        $dest->update_stats_from_disc;
        
        unless ($success) {
            unless ($d_existed) {
                $dest->remove;
            }
            $self->throw("copy of $sp => $dp failed: $error");
        }
        else {
            # check md5s match
            my $smd5 = $self->md5;
            unless ($smd5) {
                $self->update_md5;
                $smd5 = $self->md5;
            }
            $dest->update_md5;
            my $dmd5 = $dest->md5;
            unless ($dmd5 eq $smd5) {
                $dest->remove;
                $self->throw("copied $sp => $dp, but the md5s didn't match!");
            }
            
            $dest->add_metadata($self->metadata);
            
            # update the graph db if the source was in the graph
            $self->_vrpipe_schema->copy_filesystemelement($self->path->stringify, $dest->path->stringify);
            
            # if the source file was in a world-readable directory, make sure
            # the dest file's parent dirs are all world-readable
            my $dir = $self->dir;
            if ($dir ne $dest->dir) {
                my @stat = stat($dir);
                if ($stat[2] & S_IROTH) {
                    $dir = $dest->dir;
                    my %dirs;
                    $dirs{$dir} = 1;
                    my $num_parents = $dir->dir_list;
                    for (1 .. $num_parents) {
                        $dir = $dir->parent;
                        $dirs{$dir} = 1;
                    }
                    open(my $pipe, "| xargs chmod -f o+rX");
                    foreach my $dir (keys %dirs) {
                        print $pipe $dir, "\n";
                    }
                    close($pipe);
                }
            }
            
            return 1;
        }
    }
    alias cp => 'copy';
    
    method check_destination_space (Dir $dir, Int $percent_free_required = 5, Bool $throw = 1) {
        # we won't use a dir if there is less than 5% disk space left on the
        # disk it is mounted on. We obviously also won't use a dir if it doesn't
        # have enough disk space to hold our file
        my $ref = dfportable($dir);
        if (defined($ref)) {
            if ($ref->{bavail} < $self->s) {
                $self->throw("There is not enough disk space available at $dir to hold " . $self->path) if $throw;
                return 0;
            }
            
            my $total_bytes  = $ref->{blocks};
            my $bytes_free   = $ref->{bfree};
            my $percent_free = 100 / $total_bytes * $bytes_free; # ($ref->{per} is user-specific)
            
            if ($percent_free < $percent_free_required) {
                $self->throw("There is not enough disk space remaining at $dir ($percent_free\%) to safely copy or move files there.") if $throw;
                return 0;
            }
        }
        
        return 1;
    }
    
    method update_stats_from_disc (PositiveInt :$retries = 1, Bool :$keep_stat_cache = 0) {
        return unless $self->protocol eq 'file:/';
        
        my $current_s     = $self->s;
        my $current_mtime = $self->mtime;
        my $path          = $self->path;
        
        my ($new_e, $new_s, $new_mtime);
        my $trys = 0;
        while (1) {
            $self->_clear_stat_cache($path) unless $keep_stat_cache;
            $new_e     = $self->check_file_existence_on_disc($path);
            $new_s     = $self->check_file_size_on_disc($path);
            $new_mtime = $self->check_mtime_on_disc($path);
            last if $new_s != $current_s;
            last if $current_mtime ne $new_mtime;
            last if ++$trys == $retries;
            sleep 1;
        }
        $self->_clear_stat_cache($path) unless $keep_stat_cache;
        
        if (!$new_s || $current_s != $new_s || $current_mtime ne $new_mtime) {
            $self->_lines(undef);
            $self->e($new_e);
            $self->s($new_s);
            $self->mtime($new_mtime);
            if ($current_s != $new_s) {
                $self->md5(undef);
            }
            $self->update;
        }
        
        my $resolve = $self->resolve;
        if ($resolve->id != $self->id) {
            $resolve->update_stats_from_disc(retries => $retries, keep_stat_cache => $keep_stat_cache);
        }
    }
    
    method lines (Bool :$raw = 0) {
        my $s = $self->s || return 0;
        
        my $lines = $self->_lines;
        if ($raw || !$lines) {
            my $type = $raw ? 'any' : $self->type;
            my $ft = VRPipe::FileType->create($type, { file => $self->path });
            $self->disconnect if $s > 640000;
            $lines = $ft->num_lines;
            unless ($raw) {
                $self->_lines($lines);
                $self->update;
            }
        }
        
        $lines || $self->throw("Failed to find any lines in " . $self->path . ", even though it has size!");
        return $lines;
    }
    
    after _lines {
        my $resolve = $self->resolve;
        if ($resolve ne $self && defined $self->{_lines}) {
            $resolve->_lines($self->{_lines}); #*** to avoid recursion we access the self hash!! rework?...
            $resolve->update;
        }
    }
    
    method num_records {
        my $s = $self->s || return 0;
        $self->disconnect if $s > 640000;
        return $self->_filetype->num_records;
    }
    
    method num_header_lines {
        return $self->_filetype->num_header_lines;
    }
    
    method header_lines {
        return $self->_filetype->header_lines;
    }
    
    method update_md5 (Str $md5?) {
        $md5 ||= $self->file_md5($self);
        $self->md5($md5);
        $self->update;
    }
    
    after md5 {
        my $resolve = $self->resolve;
        if ($resolve ne $self && defined $self->{md5}) {
            $resolve->md5($self->{md5});       #*** to avoid recursion we access the self hash!! rework?...
            $resolve->update;
        }
    }
    
    method create_fofn (ArrayRef['VRPipe::File'] $files!) {
        my $ofh = $self->openw;
        foreach my $file (@$files) {
            my $file_path = $file->path;
            print $ofh "$file_path\n";
        }
        $ofh->close;
        $self->update_stats_from_disc(retries => 3);
        $self->throw("fofn " . $self->path . " does not contain all input files") unless ($self->lines == scalar @$files);
    }
}

1;
