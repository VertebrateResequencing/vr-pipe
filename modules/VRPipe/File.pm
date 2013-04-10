
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

class VRPipe::File extends VRPipe::Persistent {
    use Devel::GlobalDestruction;
    use MooseX::Aliases;
    use File::ReadBackwards;
    use IO::Uncompress::AnyUncompress;
    use VRPipe::FileType;
    use File::Copy;
    use Cwd qw(abs_path);
    use Filesys::DfPortable;
    
    our $bgzip_magic = [37, 213, 10, 4, 0, 0, 0, 0, 0, 377, 6, 0, 102, 103, 2, 0];
    our %file_type_map = (fastq => 'fq');
    
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
    
    has 'metadata' => (
        is      => 'rw',
        isa     => 'HashRef',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => sub { {} }
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
    
    has _opened_for_writing => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    has _opened => (
        is  => 'rw',
        isa => 'Maybe[IO::File|FileHandle]'
    );
    
    method check_file_existence_on_disc (File $path?) {
        $path ||= $self->path;         # optional so that we can call this without a db connection by supplying the path
        $self->throw("no path!") unless $path;
        my $e = -e $path;
        return $e || 0;
    }
    
    method check_file_size_on_disc (File $path?) {
        $path ||= $self->path;
        
        # we want the size of the real file, not of a symlink
        if (-l $path) {
            my $start_dir = file($path)->dir;
            while (-l $path) {
                $path = file(readlink($path));
                unless ($path->is_absolute) {
                    $path = file($start_dir, $path);
                }
                $start_dir = $path->dir;
            }
        }
        my $s = -s $path;
        
        return $s ? $s : 0;
    }
    
    method check_mtime_on_disc (File $path?) {
        $path ||= $self->path;
        
        my $st    = $path->lstat;
        my $mtime = $st ? $st->mtime : 0;
        my $dt    = DateTime->from_epoch(epoch => $mtime);
        
        return $dt;
    }
    
    method _filetype_from_extension {
        my $path = $self->path;
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
    
    method add_metadata (HashRef $meta, Bool :$replace_data = 1) {
        my $transaction = sub {
            $self->lock_row($self, 1); #*** the 1 means 'no hack' which means we rely on 'READ COMMITTED' to ensure the existing_meta we get is the most up-to-date metadata, rather than the hack that means this transaction is forced to take 1 second, which is ruinous for datasource updates
            
            my $existing_meta = $self->metadata;
            
            # incase the input $meta was the same hashref as existing_meta, we need
            # a new ref or update will do nothing
            my $new_meta = {};
            while (my ($key, $val) = each %$existing_meta) {
                $new_meta->{$key} = $val;
            }
            
            while (my ($key, $val) = each %$meta) {
                # we don't always overwrite existing values
                if ($replace_data) {
                    next unless (defined $val && "$val" ne "");
                }
                else {
                    next if exists $new_meta->{$key};
                }
                
                $new_meta->{$key} = $val;
            }
            
            $self->metadata($new_meta);
            $self->update;
        };
        $self->do_transaction($transaction, "Failed to add_metadata for file " . $self->path);
        
        my $resolve = $self->resolve;
        if ($resolve ne $self) {
            $resolve->add_metadata($meta, replace_data => $replace_data);
        }
    }
    
    method openr {
        return $self->open('<');
    }
    
    method last_line {
        my $fh = $self->open('<', backwards => 1);
        my $line = <$fh>;
        close($fh);
        return $line;
    }
    
    method openw {
        return $self->open('>');
    }
    
    method open (OpenMode $mode, Str :$permissions?, Bool :$backwards?, Int :$retry = 0) {
        my $path = $self->path;
        
        $self->throw("Only modes <, > and >> are supported") unless $mode =~ /^(?:<|>)+$/;
        
        if ($mode eq '<' && !$self->e) {
            $self->update_stats_from_disc;
            # give it a second chance...
            if (!$self->e) {
                $self->throw("File '$path' does not exist, so cannot be opened for reading");
            }
        }
        
        my $fh;
        my $type = VRPipe::FileType->create($self->type, { file => $path });
        unless (defined $backwards) {
            $backwards = $type->read_backwards;
        }
        
        # set up the open command, handling compressed files automatically
        my $open_cmd = $path;
        $self->disconnect;
        my $magic = `file -bi $path`;
        ($magic) = split(';', $magic);
        if ($magic eq 'application/octet-stream' || $path =~ /\.gz$/) {
            if ($mode eq '<') {
                if ($backwards) {
                    $self->throw("Unable to read '$path' backwards when it is compressed");
                }
                
                # if it was made with Heng Li's bgzip it will be detected as a
                # gzip file, but will fail to be decompressed properly with
                # IO::Uncompress; manually detect the magic ourselves
                if ($self->check_magic($self->path, $bgzip_magic)) {
                    $open_cmd = "gunzip -c $path |";
                }
                else {
                    $fh = IO::Uncompress::AnyUncompress->new($path->stringify, AutoClose => 1);
                }
            }
            else {
                $open_cmd = "| gzip -c $mode $path";
            }
            
            open($fh, $open_cmd) unless $fh;
        }
        else {
            if ($mode eq '<' && $backwards) {
                # we'll open it with File::ReadBackwards, but first check it can
                # opened normally
                my $ok = open(my $testfh, '<', $path);
                if ($ok) {
                    close($testfh);
                    my $rs       = $type->record_separator;
                    my @frb_args = ($path);
                    push(@frb_args, $rs) if $rs;
                    tie(*BW, 'File::ReadBackwards', @frb_args);
                    $fh = \*BW;
                }
            }
            else {
                my @args = ($mode);
                push(@args, $permissions) if $permissions;
                $fh = $path->open(@args);
            }
        }
        
        if ($fh) {
            if (index($mode, '>') == 0) {
                $self->e($self->check_file_existence_on_disc);
                $self->update;
                $self->_opened_for_writing(1);
            }
        }
        else {
            if ($retry > 59 && $mode eq '<') {
                $self->throw("Failed to open '$path' after multiple retries: $!");
            }
            else {
                # we think the file exists, so sleep a second and try again
                $self->disconnect;
                sleep(1);
                $self->e($self->check_file_existence_on_disc);
                $self->update;
                return $self->open($mode, defined $permissions ? (permissions => $permissions) : (), defined $backwards ? (backwards => $backwards) : (), retry => ++$retry);
            }
        }
        
        $self->_opened($fh);
        $self->disconnect;
        return $fh;
    }
    
    method close {
        my $fh = $self->_opened || return;
        eval { CORE::close($fh); }; #*** without the eval we get [(in cleanup) Can't use an undefined value as a symbol reference at .../File/ReadBackwards.pm line 221.] ... need to fix without eval...
        $self->_opened(undef);
        if ($self->_opened_for_writing) {
            $self->update_stats_from_disc;
            $self->_opened_for_writing(0);
        }
        return 1;
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
        
        my %stepstates = map { $_->stepstate->id => $_->stepstate } VRPipe::StepOutputFile->search({ file => $self->id });
        while (my ($ss_id, $ss) = each %stepstates) {
            $ss->pipelinesetup->log_event("File->remove() called for StepOutputFile $path, " . ($worked ? 'and it worked' : 'but it failed'), dataelement => $ss->dataelement->id, stepstate => $ss->id, record_stack => 1);
        }
        
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
            $self->_check_destination_space($dp->dir);
            
            $self->disconnect;
            $success = File::Copy::move($sp, $dp);
            
            $dest->update_stats_from_disc;
            unless ($success) {
                $self->update_stats_from_disc;
                $self->throw("move of $sp => $dp failed: $!");
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
        }
    }
    
    method update_symlink (VRPipe::File $dest) {
        # if the symlink was removed from disk by a outside of
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
 
 Title   : move (alias mv)
 Usage   : my $real_file = $obj->resolve;
 Function: If this file was created as a symlink (using VRPipe::File->symlink),
           returns the VRPipe::File corresponding to the real file. If this file
           was moved somewhere else, returns the VRPipe::File corresponding to
           the current location. If neither of these is true, returns itself.
 Returns : VRPipe::File object
 Args    : n/a

=cut
    
    method resolve {
        my $links_resolved;
        my $parent = $self->parent;
        if ($parent) {
            $links_resolved = $parent->resolve;
        }
        else {
            $links_resolved = $self;
        }
        
        my $fully_resolved;
        my $moved_to = $links_resolved->moved_to;
        if ($moved_to) {
            $fully_resolved = $moved_to->resolve;
        }
        else {
            $fully_resolved = $links_resolved;
        }
        
        return $fully_resolved;
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
        $self->_check_destination_space($dest->path->dir);
        
        $self->disconnect;
        my $success = File::Copy::copy($sp, $dp);
        $dest->update_stats_from_disc;
        
        unless ($success) {
            unless ($d_existed) {
                $dest->remove;
            }
            $self->throw("copy of $sp => $dp failed: $!");
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
            return 1;
        }
    }
    alias cp => 'copy';
    
    method _check_destination_space (Dir $dir, Int $percent_free_required = 5) {
        # we won't use a dir if there is less than 5% disk space left on the
        # disk it is mounted on. We obviously also won't use a dir if it doesn't
        # have enough disk space to hold our file
        my $ref = dfportable($dir);
        if (defined($ref)) {
            $self->throw("There is not enough disk space available at $dir to hold " . $self->path) if $ref->{bavail} < $self->s;
            
            my $total_bytes  = $ref->{blocks};
            my $bytes_free   = $ref->{bfree};
            my $percent_free = 100 - (100 / $total_bytes * $bytes_free); # ($ref->{per} is user-specific)
            $self->throw("There is not enough disk space remaining at $dir ($percent_free\%) to safely copy or move files there.") if $percent_free < $percent_free_required;
        }
    }
    
    method update_stats_from_disc (PositiveInt :$retries = 1) {
        my $current_s     = $self->s;
        my $current_mtime = $self->mtime;
        my $path          = $self->path;
        
        my ($new_e, $new_s, $new_mtime);
        my $trys = 0;
        while (1) {
            $new_e     = $self->check_file_existence_on_disc($path);
            $new_s     = $self->check_file_size_on_disc($path);
            $new_mtime = $self->check_mtime_on_disc($path);
            last if $new_s != $current_s;
            last if $current_mtime ne $new_mtime;
            last if ++$trys == $retries;
            sleep 1;
        }
        
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
        if ($resolve ne $self) {
            $resolve->update_stats_from_disc(retries => $retries);
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
        
        my $records = 0;
        my $ft = VRPipe::FileType->create($self->type, { file => $self->path });
        $self->disconnect if $s > 640000;
        $records = $ft->num_records;
        
        return $records;
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
    
    sub DEMOLISH {
        return if in_global_destruction;
        shift->close;
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
