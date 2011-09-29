=head1 NAME

VRPipe::Parser::fastq - parse fastq files

=head1 SYNOPSIS

use VRPipe::Parser;

# create object, supplying fastq file
my $pars = VRPipe::Parser->create('fastq', {file => $fastq_file});

# get the array reference that will hold the most recently requested record
my $parsed_record = $pars->parsed_record();

# loop through the fastq file, getting records
while ($pars->next_record()) {
    my $id = $parsed_record->[0];
    my $seq_string = $parsed_record->[1];
    my $qual_string = $parsed_record->[2];
}

# get a list of all the sequence ids
my @ids = $pars->sequence_ids();

# get data for a specific sequence:
my $seq_string = $pars->seq($ids[0]);
my $qual_string = $pars->quality($ids[0]);

# ask if a given id is present in the fastq file:
if ($pars->exists('xyz')) {
    # ...
}

# convert a Sanger quality string to a list of quality integers
my @qualities = $pars->qual_to_ints($qual_string);

=head1 DESCRIPTION

An indexing parser for fastq files.

If you are only doing basic parsing (only wish to call next_record and none of
the other methods), you can supply true to next_record(1) to disable indexing
and reduce memory usage and increase speed.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

use VRPipe::Base;

class VRPipe::Parser::fastq with VRPipe::ParserRole {
    use Inline C => Config => FILTERS => 'Strip_POD';
    
    our %type_to_index = (sequence => 1, quality => 2);
    
    has '_saw_last_line' => (is => 'rw',
                             isa => 'Bool',
                             default => 0);
    
    has '_seq_tells' => (is => 'rw',
                         isa => 'HashRef',
                         default => sub { {} },
                         traits => ['Hash'],
                         handles => { '_set_seq_tell' => 'set',
                                      '_get_seq_tell' => 'get',
                                      '_list_seqs' => 'keys',
                                      '_seq_exists' => 'exists'});
    
=head2 parsed_record

 Title   : parsed_record
 Usage   : my $parsed_record = $obj->parsed_record()
 Function: Get the data structure that will hold the last record requested by
           next_record()
 Returns : array ref, where the elements are:
           [0] sequence id
           [1] sequence string
           [2] quality string
 Args    : n/a

=cut

=head2 next_record

 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next record from the fastq file.
 Returns : boolean (false at end of output; check the parsed_record for the
           actual result information)
 Args    : none for normal usage, true to disable indexing (saves memory, but
           breaks most other methods)

=cut
    method next_record (Bool $no_index = 0) {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # get the next entry (4 lines)
        my $tell = tell($fh);
        my $line = <$fh>;
        unless ($line) {
            $self->_saw_last_line(1);
            return;
        }
        
        # http://maq.sourceforge.net/fastq.shtml
        # @<seqname>\n<seq>\n+[<seqname>]\n<qual>\n
        # <seqname>	:=	[A-Za-z0-9_.:-]+
        # <seq>	        :=	[A-Za-z\n\.~]+
        # <qual>	:=	[!-~\n]+
        
        # seqname
        # for speed purposes, we allow anything for the seqname up to the first
        # space
        my $filename = $self->filename;
        unless (index($line, '@') == 0) {
            $self->throw("fastq file '$filename' was bad, seqname line didn't start as expected: $line");
        }
        my ($seq_name) = split(' ', $line);
        $seq_name = substr($seq_name, 1);
        
        # seq
        # apparently sequence can be split over multiple lines, but we'll take
        # the lazy way out and assume it will always be on one line. For speed,
        # we don't validate the line.
        my $seq = <$fh>;
        $seq || $self->throw("fastq file '$filename' was truncated - no sequence for $seq_name");
        chomp($seq);
        
        # seqname repeated
        $line = <$fh>;
        $line || $self->throw("fastq file '$filename' was truncated - no + line $seq_name");
        unless (index($line, '+') == 0) {
            $self->throw("fastq file '$filename' was bad, + line didn't start as expected: $line");
        }
        
        # qual
        # for speed purposes, we don't validate the line
        my $qual = <$fh>;
        $qual || $self->throw("fastq file '$filename' was truncated - no quality line for $seq_name");
        chomp($qual);
        
        # we assume that seqnames are unique in our fastqs...
        unless ($no_index) {
            $self->_set_seq_tell($seq_name => $tell);
        }
        
        my $pr = $self->parsed_record;
        $pr->[0] = $seq_name;
        $pr->[1] = $seq;
        $pr->[2] = $qual;
        
        return 1;
    }

=head2 sequence_ids

 Title   : sequence_ids
 Usage   : my @ids = $obj->sequence_ids();
 Function: Get all the sequence ids in the fastq file. NB: only works on
           seekable input.
 Returns : list of ids
 Args    : n/a

=cut
    method sequence_ids {
        my $fh = $self->fh() || return;
        $self->_save_position || return;
        
        unless ($self->_saw_last_line) {
            while ($self->next_record) {
                next;
            }
        }
        
        $self->_restore_position;
        
        return $self->_list_seqs;
    }

=head2 seq

 Title   : seq
 Usage   : my $seq = $obj->seq($id);
 Function: Get the sequence string of a particular sequence.
 Returns : string
 Args    : string (id)

=cut
    method seq (Str $id) {
        return $self->_get_seq_data($id, 'sequence');
    }

=head2 quality

 Title   : quality
 Usage   : my $quality = $obj->quality($id);
 Function: Get the quality string of a particular sequence.
 Returns : string
 Args    : string (id)

=cut
    method quality (Str $id) {
        return $self->_get_seq_data($id, 'quality');
    }
    
    method _get_seq_data (Str $seq_name, Str $type) {
        my $fh = $self->fh() || return;
        
        $self->_parse_to_seqname($seq_name) || return;
        
        $self->_save_position || return;
        
        my $tell = $self->_get_seq_tell($seq_name);
        $self->seek($tell, 0);
        $self->next_record;
        my $index = $type_to_index{$type};
        my $result = $self->parsed_record->[$index];
        
        $self->_restore_position;
        
        return $result;
    }
    
    method _parse_to_seqname (Str $seq_name, Bool $no_warning = 0) {
        my $fh = $self->fh() || return;
        $self->_save_position || return;
        
        if (! $self->_seq_exists($seq_name)) {
            if (! $self->_saw_last_line) {
                while (! $self->_seq_exists($seq_name)) {
                    $self->next_record;
                    last if $self->_saw_last_line;
                }
            }
        }
        if (! $self->_seq_exists($seq_name)) {
            $self->warn("'$seq_name' wasn't found in this fastq file") unless $no_warning;
            return;
        }
        
        $self->_restore_position;
        
        return 1;
    }
    
=head2 exists

 Title   : exists
 Usage   : if ($obj->exists($id)) { ... }
 Function: Find out if a particular sequence name is in this fastq.
 Returns : boolean
 Args    : string (id)

=cut
    method exists (Str $seq_name) {
        my $fh = $self->fh() || return;
        
        $self->_parse_to_seqname($seq_name, 1);
        
        return $self->_seq_exists($seq_name);
    }
    
    use Inline C => <<'END_C';

=head2 qual_to_ints

 Title   : qual_to_ints
 Usage   : my @qualities = $obj->qual_to_ints($quality_string);
 Function: Convert the quality string of a fastq sequence into quality integers.
           NB: this currently only works correctly for sanger (phred) quality
           strings, as found in sam files.
 Returns : list of int
 Args    : quality string

=cut

void qual_to_ints(SV* obj, char* str) {
    Inline_Stack_Vars;
    
    Inline_Stack_Reset;
    
    char c;
    int i = 0;
    while (c = str[i++]) {
        Inline_Stack_Push(sv_2mortal(newSViv(c - 33)));
    }
    
    Inline_Stack_Done;
}

END_C
}

1;
