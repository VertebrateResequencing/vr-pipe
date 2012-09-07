
=head1 NAME

VRPipe::Parser::cat - parse VRPipe concatenate files

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying a cat file
    my $pars = VRPipe::Parser->create('cat', {file => $fa_file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # loop through the file, getting records in reverse order from the end of
    # the file
    while ($pars->next_record()) {
        my @lines = @{$parsed_record};
    }

=head1 DESCRIPTION

cat files are a special B<VRPipe> filetype for files that were created using
C<concatenate()>. That method concatenates the content of 2 text files together
and adds a marker to separate the two. The resulting chunks are supposed to be
considered in reverse order.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Parser::cat with VRPipe::ParserRole {
    has 'tail_mode' => (
        is        => 'rw',
        isa       => 'Bool',
        predicate => 'tail_mode_set'
    );

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
        
        # some sort of issue with File::ReadBackwards means we can't actually
        # do <$fh> on very large files, so have to check how big the file is
        # first and turn $fh into a tail if necessary
        my $tail_mode;
        unless ($self->tail_mode_set) {
            my $file = $self->_vrpipe_file;
            if ($file->s > 100000) {
                $file->close;
                undef($fh);
                my $path = $file->path;
                open($fh, "tail -n 1002 $path |") || $self->throw("Could not get the tail of $path"); # 1002 to potentially get the whole last record if we did a 1000 max_lines concatenate
                $self->_set_fh($fh);
                $self->tail_mode(1);
                $tail_mode = 1;
                warn "The cat file $path was too large to parse properly, can only consider the last 1000 lines\n";
            }
            else {
                $self->tail_mode(0);
            }
        }
        else {
            $tail_mode = $self->tail_mode;
        }
        
        # (we're reading the file backwards if not in tail mode)
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
            return unless ($tail_mode && @lines);
        }
        
        # fill in the parsed_record
        my $pr = $self->parsed_record;
        $#$pr = -1; # truncate the ref to no elements
        my $i = 0;
        foreach my $line ($tail_mode ? (@lines) : (reverse @lines)) {
            $pr->[$i++] = $line unless ($tail_mode && $line =~ /^-{33}VRPipe--concat-{33}/);
        }
        
        return 1;
    }
}

1;
