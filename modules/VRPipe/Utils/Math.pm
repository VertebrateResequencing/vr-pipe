
=head1 NAME

VRPipe::Utils::Math - math utility functions

=head1 SYNOPSIS
    
    use VRPipe::Utils::Math;
    
    my $math_util = VRPipe::Utils::Math->new();
    
    my $median = $math_util->histogram_median({1 => 5, 2 => 10, 3 => 5});

=head1 DESCRIPTION

General utility functions for doing math/stats stuff.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Utils::Math {
    use List::Util qw(sum);

=head2 histogram_median
 
 Title   : histogram_median
 Usage   : my $median = $obj->histogram_median($histogram_hash_ref);
 Function: Get the median bin of a histogram.
 Returns : int
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut
    
    method histogram_median (HashRef $hash) {
        # find the half-way frequency
        my $total = 0;
        foreach my $freq (values %{$hash}) {
            $total += $freq;
        }
        my $half = sprintf("%0.0f", $total / 2);
        
        # find the corresponding bin
        my $median  = 0;
        my $current = 0;
        foreach my $bin (sort { $a <=> $b } keys %{$hash}) {
            $current += $hash->{$bin};
            if ($current >= $half) {
                $median = $bin;
                last;
            }
        }
        
        return $median;
    }

=head2 histogram_quartiles
 
 Title   : histogram_quartiles
 Usage   : my $median = $obj->histogram_median($histogram_hash_ref);
 Function: Get the quartliles of a histogram.
 Returns : hash with keys q1, q2, q3
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut
    
    method histogram_quartiles (HashRef $hash) {
        my %quartiles;
        my $total_freq = sum 0, values %{$hash};
        my $freq = 0;
        foreach my $key (sort { $a <=> $b } keys %{$hash}) {
            $freq += $hash->{$key};
            
            if ($freq >= 0.75 * $total_freq) {
                $quartiles{q3} = $key;
                $quartiles{q2} = $key unless defined $quartiles{q2};
                $quartiles{q1} = $key unless defined $quartiles{q1};
                last;
            }
            elsif ($freq >= 0.5 * $total_freq and !(defined $quartiles{q2})) {
                $quartiles{q2} = $key;
                $quartiles{q1} = $key unless defined $quartiles{q1};
            }
            elsif ($freq >= 0.25 * $total_freq and !(defined $quartiles{q1})) {
                $quartiles{q1} = $key;
            }
        }
        
        return %quartiles;
    }

=head2 histogram_mean
 
 Title   : histogram_mean
 Usage   : my $mean = $obj->histogram_mean($histogram_hash_ref);
 Function: Get the mean of a histogram
 Returns : float
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut
    
    method histogram_mean (HashRef $hash) {
        my ($total, $count);
        
        while (my ($k, $v) = each(%$hash)) {
            $total += $k * $v;
            $count += $v;
        }
        
        return defined $count ? $total / $count : undef;
    }

=head2 histogram_stats
 
 Title   : histogram_stats
 Usage   : my $mean = $obj->histogram_stats($histogram_hash_ref, $mean);
 Function: Get the standard deviation of a histogram
 Returns : hash of statistics, keys are:
           mean, standard_deviation, q1, q2, q3, total (=sum of values)
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins
           Optional mean of histogram, to save calculating it.

=cut
    
    method histogram_stats (HashRef $hash, Num $mean?) {
        my %stats;
        my $sd    = 0;
        my $total = 0;
        $mean = $self->histogram_mean($hash) unless $mean;
        return undef unless defined $mean;
        $stats{mean} = $mean;
        
        # we use the formula
        # sd^2 = sum( (x_i - mean)^2 ) / n
        while (my ($num, $freq) = each %{$hash}) {
            $sd += $freq * (($mean - $num)**2);
            $total += $freq;
        }
        
        $stats{standard_deviation} = ($sd / $total)**0.5;
        $stats{total}              = $total;
        my %quartiles = $self->histogram_quartiles($hash);
        foreach (qw[q1 q2 q3]) {
            $stats{$_} = $quartiles{$_};
        }
        return %stats;
    }

=head2 compare_hash_keys
 
 Title   : compare_hash_keys
 Usage   : my $same = $obj->compare_hash_keys($hash1, $hash2);
 Function: See if two hashes have the same set of keys.
 Returns : boolean
 Args    : 2 hash refs

=cut
    
    method compare_hash_keys (HashRef $hash1, HashRef $hash2) {
        my %hash1 = %{$hash1};
        my %hash2 = %{$hash2};
        
        my $same = 1;
        if (keys %hash1 != keys %hash2) {
            $same = 0;
        }
        else {
            foreach my $key (keys %hash1) {
                last unless exists $hash2{$key};
                delete $hash2{$key};
            }
            if (%hash2) {
                $same = 0;
            }
        }
        
        return $same;
    }
}

1;
