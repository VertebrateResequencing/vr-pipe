use VRPipe::Base;

class VRPipe::Parser::cat with VRPipe::ParserRole {
=head2 parsed_record

 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, where the elements are each line of output
 Args    : n/a

=cut
    
=head2 next_record

 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next report from the cat file, starting with the last and
           working backwards.
 Returns : boolean (false at end of output; check the parsed_record for the
           actual information)
 Args    : n/a

=cut
    method next_record {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # (we're reading the file backwards)
        my ($saw_marker, @lines) = (0);
        while (my $line = $self->_readline($fh)) { #*** move the getting of $fh to inside _readline? needs to be profiled
            if ($line =~ /^-{33}VRPipe--concat-{33}/) {
                if ($saw_marker) {
                    $self->_pushback($line);
                    last;
                }
                $saw_marker = 1;
                next;
            }
            chomp($line);
            push(@lines, $line);
        }
        
        unless ($saw_marker) {
            return;
        }
        
        # fill in the parsed_record
        my $pr = $self->parsed_record;
        $#$pr = -1; # truncate the ref to no elements
        my $i = 0;
        foreach my $line (reverse @lines) {
            $pr->[$i++] = $line;
        }
        
        return 1;
    }
}

1;