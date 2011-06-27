use VRPipe::Base;

class VRPipe::File extends VRPipe::Persistent {
    use MooseX::Aliases;
    use File::ReadBackwards;
    use IO::Uncompress::AnyUncompress;
    use VRPipe::FileType;
    
    our $bgzip_magic = [37, 213, 10, 4, 0, 0, 0, 0, 0, 377, 6, 0, 102, 103, 2, 0];
    
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'path' => (is => 'rw',
                   isa => AbsoluteFile, # we can't be nice and auto convert relative to absolute because alterations made by moose during construction do not affect what gets put in the db
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   handles => [qw(slurp stat lstat)]);
    
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
                isa => IntSQL[64],
                traits => ['VRPipe::Persistent::Attributes'],
                builder => 'check_file_size_on_disc');
    
    has 'md5' => (is => 'rw',
                  isa => Varchar[64],
                  traits => ['VRPipe::Persistent::Attributes'],
                  is_nullable => 1);
    
    has _opened_for_writing => (is => 'rw',
                                isa => 'Bool',
                                default => 0);
    
    has _opened => (is => 'rw',
                    isa => 'Maybe[IO::File|FileHandle]');
    
    method check_file_existence_on_disc {
        my $e = -e $self->path;
        return $e || 0;
    }
    
    method check_file_size_on_disc {
        my $s = -s $self->path;
        return $s || 0;
    }
    
    method _filetype_from_extension {
        my $path = $self->path;
        my ($type) = $path =~ /\.(\w{2,3})$/;
        return $type || 'any';
    }
    
    __PACKAGE__->make_persistent();
    
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
    method open (OpenMode $mode, Str :$permissions?, Bool :$backwards?) {
        my $path = $self->path;
        
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
                    $fh = IO::Uncompress::AnyUncompress->new($path, AutoClose => 1);
                }
            }
            else {
                $open_cmd = "| gzip -c > $path";
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
            $self->throw("Failed to open '$path': $!");
        }
        
        $self->_opened($fh);
        return $fh;
    }
    
    method close {
        my $fh = $self->_opened || return;
        close($fh);
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
        $self->path->remove;
        $self->update_stats_from_disc;
    }
    alias unlink => 'remove';
    alias rm => 'remove';
    alias delete => 'remove';
    
    method update_stats_from_disc {
        $self->e($self->check_file_existence_on_disc);
        $self->s($self->check_file_size_on_disc);
        $self->update;
    }
    
    method DEMOLISH {
        $self->close;
    }
}

1;