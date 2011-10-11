use VRPipe::Base;

class VRPipe::File extends VRPipe::Persistent {
    use Devel::GlobalDestruction;
    use MooseX::Aliases;
    use File::ReadBackwards;
    use IO::Uncompress::AnyUncompress;
    use VRPipe::FileType;
    
    our $bgzip_magic = [37, 213, 10, 4, 0, 0, 0, 0, 0, 377, 6, 0, 102, 103, 2, 0];
    our %file_type_map = (fastq => 'fq');
    
    has 'path' => (is => 'rw',
                   isa => AbsoluteFile, # we can't be nice and auto convert relative to absolute because alterations made by moose during construction do not affect what gets put in the db
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   handles => [qw(slurp stat lstat basename dir)]);
    
    has 'type' => (is => 'rw',
                   isa => FileType,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   builder => '_filetype_from_extension');
    
    has 'e' => (is => 'rw',
                isa => 'Bool',
                traits => ['VRPipe::Persistent::Attributes'],
                builder => 'check_file_existence_on_disc');
    
    has 's' => (is => 'rw',
                isa => IntSQL[16],
                traits => ['VRPipe::Persistent::Attributes'],
                builder => 'check_file_size_on_disc');
    
    has 'md5' => (is => 'rw',
                  isa => Varchar[64],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_nullable => 1);
    
    has '_lines' => (is => 'rw',
                    isa => IntSQL[16],
                    traits => ['VRPipe::Persistent::Attributes'],
                    is_nullable => 1);
    
    has 'metadata' => (is => 'rw',
                       isa => 'HashRef',
                       traits => ['VRPipe::Persistent::Attributes'],
                       default => sub { {} });
    
    has _opened_for_writing => (is => 'rw',
                                isa => 'Bool',
                                default => 0);
    
    has _opened => (is => 'rw',
                    isa => 'Maybe[IO::File|FileHandle]');
    
    method check_file_existence_on_disc (File $path?) {
        $path ||= $self->path; # optional so that we can call this without a db connection by supplying the path
        my $e = -e $path;
        return $e || 0;
    }
    
    method check_file_size_on_disc (File $path?) {
        $path ||= $self->path;
        
        # we want the size of the real file, not of a symlink
        while (-l $path) {
            $path = readlink($path);
        }
        my $s = -s $path;
        
        return $s ? $s : 0;
    }
    
    method _filetype_from_extension {
        my $path = $self->path;
        my ($type) = $path =~ /\.([^.]*)$/;
        $type ||= 'any';
        $type = lc($type);
        if ($type eq 'gz') {
            $path =~ s/\.gz$//;
            $type = VRPipe::File->get(path => $path)->type;
        }
        if (exists $file_type_map{$type}) {
            $type = $file_type_map{$type};
        }
        eval "require VRPipe::FileType::$type;";
        return $@ ? 'any' : $type;
    }
    
    __PACKAGE__->make_persistent();
    
    method add_metadata (HashRef $meta, Bool :$replace_data = 1) {
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
        
        if ($mode eq '<' && ! $self->e) {
            $self->throw("File '$path' does not exist, so cannot be opened for reading");
        }
        
        my $fh;
        my $type = VRPipe::FileType->create($self->type, {file => $path});
        unless (defined $backwards) {
            $backwards = $type->read_backwards;
        }
        
        # set up the open command, handling compressed files automatically
        my $open_cmd = $path;
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
                my $rs = $type->record_separator;
                my @frb_args = ($path);
                push(@frb_args, $rs) if $rs;
                tie(*BW, 'File::ReadBackwards', @frb_args);
                $fh = \*BW;
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
                sleep(1);
                return $self->open($mode,
                                   defined $permissions ? (permissions => $permissions) : (),
                                   defined $backwards ? (backwards => $backwards) : (),
                                   retry => ++$retry);
            }
        }
        
        $self->_opened($fh);
        return $fh;
    }
    
    method close {
        my $fh = $self->_opened || return;
        eval {CORE::close($fh);}; #*** without the eval we get [(in cleanup) Can't use an undefined value as a symbol reference at .../File/ReadBackwards.pm line 221.] ... need to fix without eval...
        $self->_opened(undef);
        if ($self->_opened_for_writing) {
            $self->update_stats_from_disc;
            $self->_opened_for_writing(0);
        }
        return 1;
    }
    
    method touch {
        $self->path->touch;
        unless ($self->e) {
            $self->update_stats_from_disc;
        }
    }
    
    method remove {
        my $path = $self->path;
        $self->path->remove;
        $self->update_stats_from_disc;
    }
    alias unlink => 'remove';
    alias rm => 'remove';
    
    method update_stats_from_disc (PositiveInt :$retries = 1) {
        my $current_s = $self->s;
        my $path = $self->path;
        $self->disconnect;
        
        my ($new_e, $new_s);
        my $trys = 0;
        while (1) {
            $new_e = $self->check_file_existence_on_disc($path);
            $new_s = $self->check_file_size_on_disc($path);
            last if $new_s != $current_s;
            last if ++$trys == $retries;
            sleep 1;
        }
        
        if (! $new_s || $current_s != $new_s) {
            $self->_lines(undef);
            $self->e($new_e);
            $self->s($new_s);
            $self->update;
        }
    }
    
    method lines {
        my $s = $self->s || return 0;
        
        my $lines = $self->_lines;
        unless ($lines) {
            my $ft = VRPipe::FileType->create($self->type, {file => $self->path});
            $lines = $ft->num_lines;
            $self->_lines($lines);
            $self->update;
        }
        
        $lines || $self->throw("Failed to find any lines in ".$self->path.", even though it has size!");
        return $lines;
    }
    
    method num_records {
        my $s = $self->s || return 0;
        
        my $records = 0;
        my $ft = VRPipe::FileType->create($self->type, {file => $self->path});
        $records = $ft->num_records;
        
        return $records;
    }
    
    method update_md5 (Str $md5?) {
        $md5 ||= $self->file_md5($self);
        $self->md5($md5);
        $self->update;
    }
    
    sub DEMOLISH {
        return if in_global_destruction;
        shift->close;
    }
}

1;