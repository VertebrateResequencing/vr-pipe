#!/usr/bin/env perl
use strict;
use warnings;

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

use Getopt::Long;
use VRPipe::Persistent::Schema;
use VRPipe::Parser;

my $help = 0;
GetOptions("help" => \$help);

my $old_si = shift || ($help = 1);
my $new_si = shift || ($help = 1);

if ($help) {
    print <<HELP;
See what changed between 2 DCC sequence index files:
diff_sequence_index.pl old.sequence.index new.sequence.index

HELP
    exit;
}

# parse all the data out of the sequence index files
my %of_interest = (study => 3,
                   study_name => 4,
                   center_name => 5,
                   sample_id => 8,
                   sample => 9,
                   population => 10,
                   platform => 12,
                   library => 14,
                   insert_size => 17,
                   withdrawn => 20,
                   reads => 23,
                   bases => 24);
my (%old_data, %new_data);
parse($old_si, \%old_data);
parse($new_si, \%new_data);

# for lanes in old, see if there are any differences in new
my %not_in_new;
my %withdrawn_stats;
delete $of_interest{withdrawn};
my $changed = 0;
while (my ($lane, $old_data) = each %old_data) {
    my $new_data = $new_data{$lane};
    unless ($new_data) {
        $not_in_new{$lane} = 1;
        next;
    }
    
    my $was_withdrawn = $old_data->{withdrawn} ? 1 : 0;
    my $now_withdrawn = $new_data->{withdrawn} ? 1 : 0;
    if ($now_withdrawn != $was_withdrawn) {
        if ($now_withdrawn) {
            $withdrawn_stats{got_withdrawn}++;
        }
        else {
            $withdrawn_stats{got_reinstated}++;
            next;
        }
    }
    next if $now_withdrawn;
    
    my %diffs;
    my $old_md5s = join(',', sort keys %{$old_data->{md5s}});
    my $new_md5s = join(',', sort keys %{$new_data->{md5s}});
    if ($old_md5s ne $new_md5s) {
        $diffs{md5s} = [$old_md5s, $new_md5s];
    }
    foreach my $key (keys %of_interest) {
        my $old = $old_data->{$key};
        my $new = $new_data->{$key};
        if ($old ne $new) {
            $diffs{$key} = [$old, $new];
        }
    }
    
    if (keys %diffs) {
        $changed++;
        print "Lane $lane has changed:\n";
        foreach my $key (sort keys %diffs) {
            my ($old, $new) = @{$diffs{$key}};
            print "\t$key = $old => $new\n";
        }
    }
}

my $num_gone = keys %not_in_new;
if ($num_gone) {
    print "\n$num_gone lanes have disappeared completely in the new sequence.index file:\n";
    print join(' ', keys %not_in_new), "\n";
}

# see what's completely new
my %not_in_old;
my %samples_not_in_old;
while (my ($lane, $new_data) = each %new_data) {
    my $old_data = $old_data{$lane};
    unless ($old_data) {
        $not_in_old{$lane} = 1;
        $samples_not_in_old{$new_data->{sample}} = 1;
    }
}
my $num_new = keys %not_in_old;
if ($num_new) {
    my $num_new_samples = keys %samples_not_in_old;
    print "\n$num_new lanes ($num_new_samples samples) are completely new in the new sequence.index file:\n";
    print join(' ', keys %not_in_old), "\n";
}

# report on withdrawn
my $withdrawn = $withdrawn_stats{got_withdrawn} || 0;
my $reinstated = $withdrawn_stats{got_reinstated} || 0;
print "\n$withdrawn fastqs were withdrawn and $reinstated fastqs were reinstated\n";

exit;

sub parse {
    my ($si_file, $hash) = @_;
    my $pars = VRPipe::Parser->create('sequence_index', {file => $si_file});
    my $pr = $pars->parsed_record;
    
    while ($pars->next_record) {
        my $unique_name = "$pr->[2].$pr->[25]";
        $unique_name =~ s/\s/_/g;
        my $lane_data = $hash->{$unique_name} || {};
        
        $lane_data->{md5s}->{$pr->[1]} = 1;
        
        while (my ($key, $index) = each %of_interest) {
            $lane_data->{$key} = $pr->[$index];
        }
        
        $hash->{$unique_name} = $lane_data;
    }
}