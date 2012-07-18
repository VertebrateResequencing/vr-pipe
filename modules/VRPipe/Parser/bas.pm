
=head1 NAME

VRPipe::Parser::bas - parse bas files

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying bas file
    my $pars = VRPipe::Parser->create('bas', {file => $bas_file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # loop through the file, getting records
    while ($pars->next_record()) {
        my $total_bases = $result_holder->[7];
        # etc.
    }
    
    # there are a number of methods for getting stats summed for all readgroups
    # (lines) in the file, eg:
    my $total_bases = $pars->total_bases;
    my $duplicate_bases = $pars->duplicate_bases;

=head1 DESCRIPTION

A parser for bas files, which are bam statatisc files, as produced by
C<VRPipe::Steps::bam_stats->bas()>.

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

class VRPipe::Parser::bas with VRPipe::ParserRole {
=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, where the elements are:
           [0] bam_filename
           [1] bam md5 checksum
           [2] study
           [3] sample
           [4] platform
           [5] library
           [6] readgroup
           [7] #_total_bases
           [8] #_mapped_bases
           [9] #_total_reads
           [10] #_mapped_reads
           [11] #_mapped_reads_paired_in_sequencing
           [12] #_mapped_reads_properly_paired
           [13] %_of_mismatched_bases
           [14] average_quality_of_mapped_bases
           [15] mean_insert_size
           [16] insert_size_sd
           [17] median_insert_size
           [18] insert_size_median_absolute_deviation
           [19] #_duplicate_reads
           [20] #_duplicate_bases
 Args    : n/a

=cut

=head2 next_record
 
 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next line from the bas file
 Returns : boolean (false at end of output; check the parsed_record for the
           actual information)
 Args    : n/a

=cut
    
    method next_record {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # get the next line
        my $line = <$fh> || return;
        chomp($line);
        
        my @data = split(qr/\t/, $line);
        @data || return;
        
        # old bas files didn't have md5 checksum in column 2, and oldest didn't
        # have bam filename in column 1, and newest gained a 20 and 21 column
        my $pr = $self->parsed_record;
        if (@data == 19 || @data == 20 || @data == 21) {
            for my $i (0 .. $#data) {
                $pr->[$i] = $data[$i];
            }
        }
        elsif (@data == 18) {
            my $j = 0;
            for my $i (0 .. $#data) {
                if ($i == 1) {
                    $j++;
                }
                $pr->[$j] = $data[$i];
                $j++;
            }
        }
        elsif (@data == 17) {
            my $j = 2;
            for my $i (0 .. $#data) {
                $pr->[$j] = $data[$i];
                $j++;
            }
        }
        else {
            $self->throw("Unexpected number of columns (" . scalar(@data) . "); is this really a bas file?\n$line");
        }
        
        # header?
        if ($pr->[18] eq 'insert_size_median_absolute_deviation') {
            # initialise everything to unknown/0 so that if user looks at pr
            # without checking that next_record returned true, they don't
            # get the header values
            for my $i (0 .. 6) {
                $pr->[$i] = 'unknown';
            }
            for my $i (7 .. 20) {
                $pr->[$i] = 0;
            }
            $self->_set_header_parsed() unless $self->_header_parsed();
            return $self->next_record;
        }
        
        return 1;
    }

=head2 total_reads
 
 Title   : total_reads
 Usage   : my $total_reads = $obj->total_reads();
 Function: Get the total reads of all readgroups reported in the bas file.
 Returns : int
 Args    : n/a

=cut
    
    method total_reads {
        return $self->_total(9);
    }
    
    method _total (Int $index) {
        my $fh = $self->fh() || return;
        $self->_save_position || return;
        
        $self->_seek_first_record();
        
        my $total = 0;
        my $pr    = $self->parsed_record;
        while ($self->next_record) {
            $total += $pr->[$index];
        }
        
        $self->_restore_position;
        
        return $total;
    }

=head2 mapped_reads
 
 Title   : mapped_reads
 Usage   : my $mapped_reads = $obj->mapped_reads();
 Function: Get the total mapped reads of all readgroups reported in the bas
           file.
 Returns : int
 Args    : n/a

=cut
    
    method mapped_reads {
        return $self->_total(10);
    }

=head2 total_bases
 
 Title   : total_bases
 Usage   : my $total_bases = $obj->total_bases();
 Function: Get the total bases of all readgroups reported in the bas
           file.
 Returns : int
 Args    : n/a

=cut
    
    method total_bases {
        return $self->_total(7);
    }

=head2 mapped_bases
 
 Title   : mapped_bases
 Usage   : my $mapped_bases = $obj->mapped_bases();
 Function: Get the total mappedbases of all readgroups reported in the bas
           file.
 Returns : int
 Args    : n/a

=cut
    
    method mapped_bases {
        return $self->_total(8);
    }

=head2 percent_mapped
 
 Title   : percent_mapped
 Usage   : my $percent_mapped = $obj->percent_mapped();
 Function: Get the percent of mapped reads across all readgroups reported in
           the bas file
 Returns : float
 Args    : n/a

=cut
    
    method percent_mapped {
        my $total  = $self->total_reads;
        my $mapped = $self->mapped_reads;
        return (100 / $total) * $mapped;
    }

=head2 duplicate_reads
 
 Title   : duplicate_reads
 Usage   : my $duplicate_reads= $obj->duplicate_reads();
 Function: Get the total number of duplicate reads of all readgroups reported in
           the bas file.
 Returns : int
 Args    : n/a

=cut
    
    method duplicate_reads {
        return $self->_total(19);
    }

=head2 duplicate_bases
 
 Title   : duplicate_bases
 Usage   : my $duplicate_bases = $obj->duplicate_bases();
 Function: Get the total duplicate bases of all readgroups reported in the bas
           file.
 Returns : int
 Args    : n/a

=cut
    
    method duplicate_bases {
        return $self->_total(20);
    }
}

1;
