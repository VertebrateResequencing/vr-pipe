
=head1 NAME

VRPipe::Utils::GenomeChunking - genome chunking functions

=head1 SYNOPSIS
    
    use VRPipe::Utils::GenomeChunking;
    
    my $chunk_util = VRPipe::Utils::GenomeChunking->new();
    
    my $chunks = $chunk_util->chunks({1 => 5, 2 => 10, 3 => 5});

=head1 DESCRIPTION

General utility functions for doing math/stats stuff.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Utils::GenomeChunking {

=head2 chunks
 
 Title   : chunks
 Usage   : my $chunks = $obj->chunks(reference_index => $reference_index, chunk_size => 100000);
 Function: Splits a reference .
 Returns : ArrayRef of chunk HashRefs. Chunks HashRef consists of keys chrom, from and to and seq_no. Optionally, male_ploidy and female_ploidy.
 Args    : 

=cut
    
    method chunks (Str|File :$reference_index!, Int :$chunk_size! = 1000000, Int :$chunk_overlap! = 0, ArrayRef :$chroms = [], HashRef :$ploidy = {}) {
        my $ploidy_default = $$ploidy{default} || 2;
        
        my $pars = VRPipe::Parser->create('fai', { file => $reference_index });
        my $parsed_record = $pars->parsed_record();
        
        my @seqs = @$chroms;
        unless (@seqs) {
            @seqs = $pars->seq_ids;
        }
        
        my @chr_regions;
        while ($pars->next_record()) {
            my ($chr, $from, $to);
            for my $regex (@seqs) {
                my $seq = $parsed_record->{SN};
                if (!($seq =~ /^([cC]hr|[cC]hrom)?($regex)$/i)) { next; }
                $chr  = $seq;
                $from = 1;
                $to   = $parsed_record->{LN};
                last;
            }
            if (!defined $chr) { next; }
            
            if (!exists($$ploidy{$chr})) {
                push @chr_regions, { chrom => $chr, from => $from, to => $to, female => $ploidy_default, male => $ploidy_default };
                next;
            }
            
            # Split the chunks as necessary
            if (exists $$ploidy{$chr} && ref($$ploidy{$chr}) ne 'ARRAY') { $self->throw("Ploidy definition not well defined for $chr regions.\n"); }
            foreach my $reg (sort { $$a{from} <=> $$b{to} } @{ $$ploidy{$chr} }) {
                my $start  = $$reg{from};
                my $end    = $$reg{to};
                my $female = defined $$reg{F} ? $$reg{F} : $ploidy_default;
                my $male   = defined $$reg{M} ? $$reg{M} : $ploidy_default;
                if ($start > $from) {
                    push @chr_regions, { chrom => $chr, from => $from, to => $start - 1, female => $ploidy_default, male => $ploidy_default };
                }
                push @chr_regions, { chrom => $chr, from => $start, to => $end, female => $female, male => $male };
                $from = $end + 1;
            }
            
            if ($from < $to) {
                push @chr_regions, { chrom => $chr, from => $from, to => $to, female => $ploidy_default, male => $ploidy_default };
            }
        }
        
        my @chunks;
        my $seq_no = 1;
        for my $region (@chr_regions) {
            my $pos     = $$region{from};
            my $end_pos = $$region{to};
            my $female  = $$region{female};
            my $male    = $$region{male};
            while ($pos < $end_pos) {
                my $from = $pos;
                my $to   = $from + $chunk_size - 1;
                
                if ($to > $end_pos) { $to = $end_pos; }
                
                if (keys %$ploidy) {
                    push @chunks, { chrom => $$region{chrom}, from => $from, to => $to, female_ploidy => $female, male_ploidy => $male, seq_no => $seq_no };
                }
                else {
                    push @chunks, { chrom => $$region{chrom}, from => $from, to => $to, seq_no => $seq_no };
                }
                ++$seq_no;
                
                $pos += $chunk_size - $chunk_overlap;
                if ($pos < 1) { $self->throw("The split size too small [$chunk_size]?\n"); }
            }
        }
        return \@chunks;
    }
}

1;
