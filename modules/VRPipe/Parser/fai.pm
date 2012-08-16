
=head1 NAME

VRPipe::Parser::fai - parse fasta index files

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying fasta index file
    my $pars = VRPipe::Parser->create('fai', {file => $fai_file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # loop through the fai file, getting records
    while ($pars->next_record()) {
        # check $parsed_record for desired info, eg:
        my $length = $parsed_record->{LN}
    }
    
    # get a hash of id => length
    my %seq_lengths = $pars->seq_lengths;
    
    # get a list of sequence ids
    my @ids = $pars->seq_ids;
    
    # get the total length of all sequences
    my $total_basespairs = $pars->total_length;

=head1 DESCRIPTION

A parser for fasta index (.fai) files.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Parser::fai with VRPipe::ParserRole {
    has '_saw_last_line' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    has 'parsed_record' => (
        is      => 'ro',
        isa     => 'HashRef',
        default => sub { {} }
    );

=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record = $obj->parsed_record()
 Function: Get the data structure that will hold the last record requested by
           next_record()
 Returns : hash ref, where the keys are:
           SN (sequence identifier)
           LN (sequence length)
 Args    : n/a

=cut

=head2 next_record
 
 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next record from the fai file.
 Returns : boolean (false at end of output; check the parsed_record for the
           actual result information)
 Args    : n/a

=cut
    
    method next_record () {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # get the next line
        my $line = <$fh> || return;
        
        # 1       249250621       52      60      61
        # 2       243199373       253404903       60      61
        # 3       198022430       500657651       60      61
        
        chomp($line);
        my @data = split(qr/\t/, $line);
        @data || return;
        
        my $pr = $self->parsed_record;
        
        foreach my $key (keys %{$pr}) {
            delete $pr->{$key};
        }
        
        $pr->{SN} = $data[0];
        $pr->{LN} = $data[1];
        
        return 1;
    }

=head2 seq_lengths
 
 Title   : seq_lengths
 Usage   : my %seq_lengths = $obj->seq_lengths;
 Function: Get all the sequence lengths.
 Returns : hash (keys as sequence ids, values as int lengths)
 Args    : n/a

=cut
    
    method seq_lengths {
        my %seq_lengths;
        
        $self->_save_position;
        $self->_seek_first_record;
        
        my $pr = $self->parsed_record;
        while ($self->next_record) {
            $seq_lengths{ $pr->{SN} } = $pr->{LN};
        }
        
        $self->_restore_position();
        return %seq_lengths;
    }

=head2 seq_ids
 
 Title   : seq_ids
 Usage   : my @seq_ids = $obj->seq_ids;
 Function: Get all the sequence ids.
 Returns : list of strings
 Args    : n/a

=cut
    
    method seq_ids {
        my @seq_ids;
        
        $self->_save_position;
        $self->_seek_first_record;
        
        my $pr = $self->parsed_record;
        while ($self->next_record) {
            push(@seq_ids, $pr->{SN});
        }
        
        $self->_restore_position();
        return @seq_ids;
    }

=head2 total_length
 
 Title   : total_length
 Usage   : my $total_length = $obj->total_length;
 Function: Get the total length of all sequences in the fai file
 Returns : int (bp)
 Args    : n/a

=cut
    
    method total_length {
        $self->_save_position;
        $self->_seek_first_record;
        
        my $pr     = $self->parsed_record;
        my $length = 0;
        while ($self->next_record) {
            $length += $pr->{LN} || 0;
        }
        
        $self->_restore_position();
        
        return $length;
    }
}

1;
