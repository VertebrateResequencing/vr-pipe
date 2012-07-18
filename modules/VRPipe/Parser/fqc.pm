
=head1 NAME

VRPipe::Parser::fqc - parse fastqcheck files

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying a fastqcheck file
    my $pars = VRPipe::Parser->create('fqc', {file => $file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # loop through the file, getting records
    while ($pars->next_record()) {
        my $position = $parsed_record->[0];
        # etc.
    }
    
    # get the summary stats of interest:
    my $reads = $pars->num_sequences();
    my $bps = $pars->total_length();
    # etc.

=head1 DESCRIPTION

fastqcheck files are produced by the B<fastqcheck> executable and give
statistics about fastq files.

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

class VRPipe::Parser::fqc with VRPipe::ParserRole {
    has 'num_sequences' => (is     => 'ro',
                            isa    => PositiveInt,
                            coerce => 1,
                            writer => '_num_sequences');
    
    has 'total_length' => (is     => 'ro',
                           isa    => PositiveInt,
                           coerce => 1,
                           writer => '_total_length');
    
    has 'avg_length' => (is     => 'ro',
                         isa    => 'Num',
                         writer => '_avg_length');
    
    has 'max_length' => (is     => 'ro',
                         isa    => PositiveInt,
                         coerce => 1,
                         writer => '_max_length');
    
    has 'standard_deviations' => (is     => 'ro',
                                  isa    => 'ArrayRef',
                                  writer => '_standard_deviations');

=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, where the elements are:
           [0] base position (starting at position 0 == 'total')
           [1] A %
           [2] C %
           [3] G %
           [4] T %
           [5] N %
           [6] count of bases at position [0] with quality 0
           [7] count of bases at position [0] with quality 1
           ... etc. for quality up to 40 (or higher?)
           [-1] the 'AQ' which is NOT the average quality... but what is it?!
 Args    : n/a

=cut

=head2 next_record
 
 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next line from the fastqcheck file.
 Returns : boolean (false at end of output; check the parsed_record for the
           actual information)
 Args    : n/a

=cut
    
    method next_record {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # get the next line
        my $line = <$fh> || return;
        
        my @data = split(qr/\s+/, $line);
        @data || return;
        
        my $pr = $self->parsed_record;
        
        my $r = 0;
        for my $i (0 .. $#data) {
            my $datum = $data[$i];
            if ($i == 0) {
                next if $datum eq 'base';
                $datum = 0 if $datum eq 'Total';
            }
            $pr->[$r++] = $datum;
        }
        
        return 1;
    }

=head2 num_sequences
 
 Title   : num_sequences
 Usage   : my $num_of_sequences = $obj->num_sequences();
 Function: Get the number of sequences that the fastqcheck file is summarising.
 Returns : int
 Args    : n/a

=cut

=head2 total_length
 
 Title   : total_length
 Usage   : my $total_length_of_all_sequences = $obj->total_length();
 Function: Get the total length of sequences that the fastqcheck file is
           summarising.
 Returns : int
 Args    : n/a

=cut

=head2 avg_length
 
 Title   : avg_length
 Usage   : my $average_length_of_a_sequence = $obj->avg_length();
 Function: Get the average length of a sequence.
 Returns : number
 Args    : n/a

=cut

=head2 max_length
 
 Title   : max_length
 Usage   : my $length_of_longest_sequence = $obj->max_length();
 Function: Get the length of the longest sequence.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

=head2 standard_deviations
 
 Title   : standard_deviations
 Usage   : my ($total_sd, $per_base_sd) = @{$obj->standard_deviations()};
 Function: Get the standard deviations at 0.25; total and per-base.
 Returns : arrayref of two percentages
 Args    : n/a

=cut
    
    method _get_header {
        my $fh = $self->fh() || return;
        return 1 if $self->_header_parsed();
        
        my $saw = 0;
        while (<$fh>) {
            if (/^(\d+) sequences, (\d+) total length, (\d+\.\d+) average, (\d+) max/) {
                $self->_num_sequences($1);
                $self->_total_length($2);
                $self->_avg_length($3);
                $self->_max_length($4);
                $saw++;
            }
            elsif (/^Standard deviations at 0.25:\s+total\s+(\d+\.\d+) %, per base\s+(\d+\.\d+) %/) {
                $self->_standard_deviations([$1, $2]);
                $saw++;
            }
            elsif (/^\s+A\s+C/) {
                $saw++;
                last;
            }
        }
        
        if ($saw == 3) {
            $self->_set_header_parsed();
            return 1;
        }
        
        $self->throw("Unable to parse header before first result - is this a fastqcheck file?");
    }

=head2 avg_base_quals
 
 Title   : avg_base_quals
 Usage   : my ($bases, $quals) = $pars->avg_base_quals();
 Function: Get the average qualities for each base.
 Returns : list with two array references
 Args    : n/a

=cut
    
    method avg_base_quals {
        $self->_save_position || return;
        $self->_seek_first_record();
        my $rh = $self->parsed_record();
        
        $self->next_record();
        
        my (@bases, @quals);
        while ($self->next_record()) {
            push @bases, $rh->[0];
            push @quals, $self->_avg_base_qual;
        }
        
        $self->_restore_position;
        
        return (\@bases, \@quals);
    }

=head2 avg_qual
 
 Title   : avg_qual
 Usage   : my $avg_qual = $pars->avg_qual();
 Function: Get the average quality of all bases.
 Returns : float
 Args    : n/a

=cut
    
    method avg_qual {
        $self->_save_position || return;
        $self->_seek_first_record();
        
        $self->next_record();
        
        my $avg_qual = $self->_avg_base_qual;
        
        $self->_restore_position;
        
        return $avg_qual;
    }
    
    method _avg_base_qual {
        my $rh = $self->parsed_record();
        
        my $sum    = 0;
        my $nvals  = 0;
        my $nquals = $#{$rh} - 1;
        for my $i (6 .. $nquals) {
            $nvals += $rh->[$i];
            $sum += $rh->[$i] * ($i - 6);
        }
        my $avg_qual = $nvals ? $sum / $nvals : 0;
        
        return $avg_qual;
    }
}

1;
