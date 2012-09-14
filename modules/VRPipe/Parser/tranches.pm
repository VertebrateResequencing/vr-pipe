
=head1 NAME

VRPipe::Parser::tranches - parse tranches file produced by GATK VQSR

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying bas file
    my $pars = VRPipe::Parser->create('tranches', {file => $tranches_file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # e.g. loop through the file, to find the minimum VQSLOD score corresponding to a
    # given truth sensitivity cutoff
    
    my $ts_filter_level = 99.85;
    my $min_vqslod = 0;
    while ($pars->next_record()) {
        next unless $parsed_record->{targetTruthSensitivity} < $ts_filter_level;
        $min_vqslod = $min_vqslod > $parsed_record->{minVQSLod} ? $min_vqslod : $parsed_record->{minVQSLod});
    }
    
=head1 DESCRIPTION

A parser for tranches files, as produced by the GATK Variant Quality Score 
Recalibration method.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Parser::tranches with VRPipe::ParserRole {

    has 'parsed_record' => (
        is      => 'ro',
        isa     => 'HashRef',
        default => sub { {} }
    );

=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : hash ref, where the keys are:
           - targetTruthSensitivity
           - numKnown
           - numNovel
           - knownTiTv
           - novelTiTv
           - minVQSLod
           - filterName
           - accessibleTruthSites
           - callsAtTruthSites
           - truthSensitivity
           - mode (if Version 5)
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
        
        # ignore header lines
        while (index($line, '#') == 0) {
            $line = <$fh> || return;
        }
        
        chomp($line);
        
        my @data = split(qr/,/, $line);
        @data || return;
        
        # # Variant quality score tranches file
        # # Version number 4
        # targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,accessibleTruthSites,callsAtTruthSites,truthSensitivity
        # 95.00,12709365,3437394,2.2470,2.1520,3.8133,TruthSensitivityTranche0.00to95.00,1422549,1351421,0.9500
        # 95.10,12773707,3576561,2.2473,2.1511,3.6949,TruthSensitivityTranche95.00to95.10,1422549,1352844,0.9510
        # 95.20,12837350,3704149,2.2477,2.1506,3.5836,TruthSensitivityTranche95.10to95.20,1422549,1354266,0.9520
        
        # # Variant quality score tranches file
        # # Version number 5
        # targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity
        # 90.00,168,0,2.5745,0.0000,2.4829,VQSRTrancheSNP0.00to90.00,SNP,165,148,0.8970
        # 99.00,963,0,2.4149,0.0000,1.6421,VQSRTrancheSNP90.00to99.00,SNP,165,163,0.9879
        # 99.90,1245,1111,2.3378,2.6426,0.7663,VQSRTrancheSNP99.00to99.90,SNP,165,164,0.9939
        # 100.00,1315,5056,2.2630,1.9935,-486.3780,VQSRTrancheSNP99.90to100.00,SNP,165,165,1.0000

        my $pr = $self->parsed_record;
        
        foreach my $key (keys %{$pr}) {
            delete $pr->{$key};
        }

        my @columns;
        if (@data == 10) {
            @columns = (qw(targetTruthSensitivity numKnown numNovel knownTiTv novelTiTv minVQSLod filterName accessibleTruthSites callsAtTruthSites truthSensitivity));
            foreach my $name (@columns) {
                $pr->{$name} = shift @data;
            }
        }
        elsif (@data == 11) {
            @columns = (qw(targetTruthSensitivity numKnown numNovel knownTiTv novelTiTv minVQSLod filterName model accessibleTruthSites callsAtTruthSites truthSensitivity));
            foreach my $name (@columns) {
                $pr->{$name} = shift @data;
            }
        }
        else {
            $self->throw("Unexpected number of columns (" . scalar(@data) . "); is this really a tranches file?\n$line");
        }
        
        # header?
        if ($pr->{targetTruthSensitivity} eq 'targetTruthSensitivity') {
            # initialise everything to unknown/0 so that if user looks at pr
            # without checking that next_record returned true, they don't
            # get the header values
            foreach my $name (@columns) {
                $pr->{$name} = 0;
            }
            $pr->{targetTruthSensitivity} = 'unknown';
            $pr->{model} = 'unknown' if exists $pr->{model};
            $self->_set_header_parsed() unless $self->_header_parsed();
            return $self->next_record;
        }
        
        return 1;
    }

}

1;
