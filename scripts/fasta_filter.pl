#!/usr/bin/env perl
use strict;
use warnings;

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