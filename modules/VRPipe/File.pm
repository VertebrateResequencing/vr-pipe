use VRPipe::Base;

class VRPipe::File extends VRPipe::Persistent {
    use MooseX::Aliases;
    
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
                   traits => ['VRPipe::Persistent::Attributes']);
    
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
                    isa => 'Maybe[IO::File]');
    
    method check_file_existence_on_disc {
        my $e = -e $self->path;
        return $e || 0;
    }
    
    method check_file_size_on_disc {
        my $s = -s $self->path;
        return $s || 0;
    }
    
    __PACKAGE__->make_persistent();
    
    method openr {
        return $self->open('r');
    }
    method openw {
        return $self->open('w');
    }
    method open (Str $mode, Str $permissions?) {
        my $path = $self->path;
        my $io_file = $path->open($mode, $permissions);
        if ($io_file) {
            if ($mode eq 'w') {
                $self->e($self->check_file_existence_on_disc);
                $self->update;
                $self->_opened_for_writing(1);
            }
        }
        else {
            $self->throw("Failed to open '$path': $!");
        }
        $self->_opened($io_file);
        return $io_file;
    }
    
    method close {
        my $io_file = $self->_opened || return;
        $io_file->close;
        $self->_opened(undef);
        if ($self->_opened_for_writing) {
            $self->s($self->check_file_size_on_disc);
            $self->update;
            $self->_opened_for_writing(0);
        }
        return 1;
    }
    
    method touch {
        $self->path->touch;
        unless ($self->e) {
            $self->e($self->check_file_existence_on_disc);
            $self->update;
        }
    }
    
    method remove {
        $self->path->remove;
        $self->e($self->check_file_existence_on_disc);
        $self->update;
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