use VRPipe::Base;

role VRPipe::ParserRole {
    has 'file' => (is => 'ro',
                   isa => VRPFileOrHandle,
                   coerce => 1,
                   required => 1);
    
    has 'type' => (is => 'ro',
                   isa => FileType,
                   lazy => 1,
                   builder => '_build_type');
    
    has 'parsed_record' => (is => 'ro',
                            isa => 'ArrayRef',
                            default => sub { [] });
    
    has 'fh' => (is => 'ro',
                 isa => AnyFileHandle,
                 lazy => 1,
                 builder => '_get_fh');
    
    has '_buffer_store' => (is => 'ro',
                            traits => ['Array'],
                            isa => 'ArrayRef[Str]',
                            default => sub { [] },
                            handles => { '_pushback' => 'push',
                                         '_getback' => 'shift',
                                         '_empty_buffer' => 'clear' });
    
    has '_vrpipe_file' => (is => 'ro',
                           isa => 'VRPipe::File',
                           writer => '_set_vrpipe_file');
    
    has '_tell' => (is => 'rw',
                    isa => 'Int');
    
    has '_current_results' => (is => 'rw',
                               isa => 'ArrayRef|HashRef');
    
    has '_header_parsed' => (is => 'ro',
                             isa => 'Int',
                             writer => '_set_header_parsed',
                             predicate => '_got_header');
    
    method _build_type {
        my $class = ref($self);
        my ($type) = $class =~ /.+::(.+)/;
        return $type;
    }
    
    method _get_fh {
        my $file = $self->file;
        if (ref($file) && (File->check($file) || $file->isa('VRPipe::File'))) {
            eval "require VRPipe::File;";
            my $vrpf = $file->isa('VRPipe::File') ? $file : VRPipe::File->get(path => $file->absolute, type => $self->type);
            $self->_set_vrpipe_file($vrpf); # because on destruction, $vrpf will close the opened filehandle
            return $vrpf->openr();
        }
        else {
            return $file;
        }
    }
    
    requires 'next_record';
    
    method _readline (AnyFileHandle $fh) {
        my $line = $self->_getback;
        unless (defined $line) {
            $line = <$fh>;
        }
        return $line;
    }
    
=head2 _save_position

 Title   : _save_position
 Usage   : $self->_save_position()
 Function: Internal method for parser authors. Saves the current filehandle
           position; for use before seeking.
 Returns : boolean (true if this is a seekable filehandle)
 Args    : n/a

=cut
    method _save_position {
        my $fh = $self->fh() || return;
        
        my $tell = tell($fh);
        if ($tell == -1) {
            $self->warn("this parsing method doesn't work on piped input");
            return;
        }
        my $current_results = ref($self->parsed_record) eq 'ARRAY' ? [@{$self->parsed_record}] : {%{$self->parsed_record}};
        
        $self->_tell($tell);
        $self->_current_results($current_results);
        
        return 1;
    }
    
=head2 _get_header

 Title   : _get_header
 Usage   : $self->_get_header()
 Function: Internal method for parser authors. This does nothing. If your file
           format has a header, you should implement it by overriding this.
           In your implementation, _set_header_parsed() should be called after
           successfully parsing your header.
           It should return true if the header was parsed, false if not.
 Returns : boolean
 Args    : n/a

=cut
    method _get_header {
        unless ($self->_header_parsed) {
            $self->_set_header_parsed(0);
        }
        return 0;
    }
    
=head2 _set_header_parsed

 Title   : _set_header_parsed
 Usage   : $self->_set_header_parsed()
 Function: Internal method for parser authors. Ask if the header of the current
           file has been parsed.
 Returns : boolean
 Args    : n/a (optionally, supply a filehandle position as from tell() that
           denotes the start of the first record after the header; by default
           it assumes the current filehandle position matches this criteria)

=cut
    around _set_header_parsed (Int $tell?) {
        my $fh = $self->fh() || return;
        unless (defined $tell) {
            $tell = tell($fh);
        }
        $self->$orig($tell);
    }
    
=head2 _seek_first_result

 Title   : _seek_first_result
 Usage   : $self->_seek_first_result()
 Function: Internal method for parser authors. Seeks back to before the first
           result (ie. after a header if present) so that next_result() will
           behave as if it was called for the first time on a file.
           You should call _save_position() first, then call this, then do your
           work, then call _restore_position().
 Returns : n/a
 Args    : n/a

=cut
    method _seek_first_result {
        my $tell;
        if ($self->_got_header) {
            $tell = $self->_header_parsed;
            if ($tell == -1) {
                $self->warn("this parsing method doesn't work on piped input");
                return;
            }
        }
        else {
            $tell = 0;
        }
        
        seek($self->fh, $tell, 0);
    }
    
=head2 _restore_position

 Title   : _restore_position
 Usage   : $self->_restore_position()
 Function: Internal method for parser authors. Restores the current filehandle
           position to the location when _save_position() was last called; for
           use after you've seeked somewhere and done some work.
           If _save_position() was called at the true start of the file, this
           actually calls _seek_first_result().
 Returns : n/a
 Args    : n/a

=cut
    method _restore_position {
        my $fh = $self->fh() || return;
        
        if ($self->_tell == 0) {
            # we might have saved position before parsing the header, but now have
            # the flag set that we've parsed the header; don't seek back before
            # the header!
            $self->_seek_first_result();
            return;
        }
        
        seek($fh, $self->{_tell}, 0);
        if (ref($self->parsed_record) eq 'ARRAY') {
            my @current_results = @{$self->_current_results};
            for my $i (0..$#current_results) {
                $self->parsed_record->[$i] = $current_results[$i];
            }
        }
        else {
            while (my ($key, $val) = each %{$self->_current_results}) {
                $self->parsed_record->{$key} = $val;
            }
        }
        
        return 1;
    }
}

1;