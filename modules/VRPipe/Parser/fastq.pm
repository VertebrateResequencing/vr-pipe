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

# convert a Sanger quality string to a list of quality integers
my @qualities = $pars->qual_to_ints($qual_string);

=head1 DESCRIPTION

A basic parser for fastq files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

use VRPipe::Base;

class VRPipe::Parser::fastq with VRPipe::ParserRole {
    use Inline C => Config => FILTERS => 'Strip_POD';
    
    has '_saw_last_line' => (is => 'rw',
                             isa => 'Bool',
                             default => 0);
    
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
 Args    : n/a

=cut
    method next_record () {
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
        
        my $pr = $self->parsed_record;
        $pr->[0] = $seq_name;
        $pr->[1] = $seq;
        $pr->[2] = $qual;
        
        return 1;
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
