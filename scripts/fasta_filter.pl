#!/usr/bin/env perl
use strict;
use warnings;

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

use Getopt::Long;
use VRPipe::Persistent::Schema;
use VRPipe::Parser;

my $help = 0;
my $retries = 3;
my ($min_length);
GetOptions("help" => \$help,
           "min_length=i" => \$min_length);

my $in_fasta = shift || ($help = 1);

if ($help) {
    print <<HELP;
Filter reads from a fasta file, producing a new (smaller) fasta file:
fasta_filter.pl [options] input.fasta > output.fasta

Options:
    --min_length <int> Minimum length of sequences that are output, in bp
HELP
    exit;
}

# create fasta parser
my $pars = VRPipe::Parser->create('fasta', {file => $in_fasta});
my $pr = $pars->parsed_record;

# loop through all sequences
my ($in_count, $out_count) = (0, 0);
while ($pars->next_record) {
    $in_count++;
    my $id = $pr->[0];
    my $seq = $pr->[1];
    
    # filter based on criteria (there's only 1 right now...)
    my $output = 0;
    if ($min_length) {
        if (length($seq) >= $min_length) {
            $output = 1;
        }
    }
    
    if ($output) {
        my @seq_lines = unpack("(A60)*", $seq);
        print '>', $id, "\n", join("\n", @seq_lines), "\n";
        $out_count++;
    }
}

warn "$in_count sequences read; $out_count sequences output\n";

exit;