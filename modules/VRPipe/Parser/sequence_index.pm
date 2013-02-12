
=head1 NAME

VRPipe::Parser::sequence_index - parse sequence.index files from the DCC

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying sequence.index file or filehandle
    my $pars = VRPipe::Parser->create('sequence_index', {file => $si_file});
    
    # get the array reference that will hold the most recently requested record
    my $parsed_record = $pars->parsed_record();
    
    # loop through the file, getting records
    while ($pars->next_record()) {
        my $fastq_file = $parsed_record->[0];
        # etc.
    }
    
    # or get specific information on a particular run id (lane)
    my $individual = $pars->lane_info('SRR005864', 'sample_name');
    
    # get all the lanes for a particular individual
    my @lanes = $pars->get_lanes('NA11994');

=head1 DESCRIPTION

sequence.index files contain metadata describing fastq files, used for the 1000
genomes project and described here: L<http://www.1000genomes.org/formats>.

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

class VRPipe::Parser::sequence_index with VRPipe::ParserRole {
    our %field_to_index = (
        FASTQ_FILE          => 0,
        MD5                 => 1,
        RUN_ID              => 2,
        STUDY_ID            => 3,
        STUDY_NAME          => 4,
        CENTER_NAME         => 5,
        SUBMISSION_ID       => 6,
        SUBMISSION_DATE     => 7,
        SAMPLE_ID           => 8,
        SAMPLE_NAME         => 9,
        POPULATION          => 10,
        EXPERIMENT_ID       => 11,
        INSTRUMENT_PLATFORM => 12,
        INSTRUMENT_MODEL    => 13,
        LIBRARY_NAME        => 14,
        RUN_NAME            => 15,
        RUN_BLOCK_NAME      => 16,
        INSERT_SIZE         => 17,
        LIBRARY_LAYOUT      => 18,
        PAIRED_FASTQ        => 19,
        WITHDRAWN           => 20,
        WITHDRAWN_DATE      => 21,
        COMMENT             => 22,
        READ_COUNT          => 23,
        BASE_COUNT          => 24,
        ANALYSIS_GROUP      => 25
    );
    
    has '_saw_last_line' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    has '_lane_tells' => (
        is      => 'rw',
        isa     => 'HashRef',
        default => sub { {} }
    );

=head2 parsed_record
 
 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, where the elements are:
           [0]  FASTQ_FILE
           [1]  MD5
           [2]  RUN_ID
           [3]  STUDY_ID
           [4]  STUDY_NAME
           [5]  CENTER_NAME
           [6]  SUBMISSION_ID
           [7]  SUBMISSION_DATE
           [8]  SAMPLE_ID
           [9]  SAMPLE_NAME
           [10] POPULATION
           [11] EXPERIMENT_ID
           [12] INSTRUMENT_PLATFORM
           [13] INSTRUMENT_MODEL
           [14] LIBRARY_NAME
           [15] RUN_NAME
           [16] RUN_BLOCK_NAME
           [17] INSERT_SIZE
           [18] LIBRARY_LAYOUT
           [19] PAIRED_FASTQ
           [20] WITHDRAWN
           [21] WITHDRAWN_DATE
           [22] COMMENT
           [23] READ_COUNT
           [24] BASE_COUNT
           [25] ANALYSIS_GROUP
 Args    : n/a

=cut

=head2 next_record
 
 Title   : next_record
 Usage   : while ($obj->next_record()) { # look in parsed_record }
 Function: Parse the next line from the sequence.index file
 Returns : boolean (false at end of output; check the parsed_record for the
           actual information)
 Args    : n/a

=cut
    
    method next_record {
        # just return if no file set
        my $fh = $self->fh() || return;
        
        # get the next line
        my $tell = tell($fh);
        my $line = <$fh>;
        unless ($line) {
            $self->_saw_last_line(1);
            return;
        }
        
        my @data = split(/\t/, $line);
        chomp($data[-1]);
        if (@data != 26 && @data != 25) {
            $self->throw("Expected 26 or 25 columns, got " . scalar @data . ": " . join(",", @data) . "\n");
            return;
        }
        
        if ($data[2] eq 'RUN_ID') {
            $self->throw("got RUN_ID!");
        }
        
        $self->_lane_tells->{ $data[2] }->{$tell} = 1;
        
        my $pr = $self->parsed_record;
        for my $i (0 .. $#data) {
            $pr->[$i] = $data[$i];
        }
        
        return 1;
    }

=head2 lane_info
 
 Title   : lane_info
 Usage   : my $sample_name = $obj->lane_info('SRR005864', 'sample_name');
 Function: Get a particular bit of info for a particular lane (run_id).
 Returns : in scalar context, the most common answer; in list context, all
           the different answers
 Args    : lane/run_id, field name (see parsed_record() for the list of valid
           field names)

=cut
    
    method lane_info (Str $lane, Str $field) {
        $field = uc($field);
        my $index = $field_to_index{$field};
        $self->throw("'$field' wasn't a valid field name") unless defined $index;
        
        my $fh = $self->fh() || return;
        $self->_save_position || return;
        
        # since a lane can appear multiple times in the file, we have to just
        # parse the whole file
        unless ($self->_saw_last_line) {
            while ($self->next_record) {
                next;
            }
        }
        my $lane_tells = $self->_lane_tells;
        if (!defined $lane_tells->{$lane}) {
            $self->warn("'$lane' wasn't a run_id in your sequence.index file");
            return;
        }
        
        my %answers;
        my $pr = $self->parsed_record;
        foreach my $tell (keys %{ $lane_tells->{$lane} }) {
            seek($fh, $tell, 0);
            $self->next_record;
            $answers{ $pr->[$index] }++;
        }
        
        $self->_restore_position;
        
        if (wantarray) {
            my @answers = sort keys %answers;
            return @answers;
        }
        else {
            my $most_common_answer;
            my $highest_count = 0;
            while (my ($answer, $count) = each %answers) {
                if ($count > $highest_count) {
                    $highest_count      = $count;
                    $most_common_answer = $answer;
                }
            }
            return $most_common_answer;
        }
    }

=head2 get_lanes
 
 Title   : get_lanes
 Usage   : my @lanes = $obj->get_lanes(sample_name => 'NA11994');
 Function: Get all the lane ids with a given property.
 Returns : string
 Args    : none for all lanes, or a single key => value pair to get lanes under
           that key, eg. sample_name => 'NA11994' to get all lanes of that
           sample.
           
           optionally, you can add one or more key value pairs where keys are
           valid field names (see result_holder() for the list) prefixed with
           'ignore_' and values are what should should be ignored (a regex).
           Eg. to ignore all withdrawn lanes:
           ignore_withdrawn => 1

=cut
    
    method get_lanes (%args) {
        $self->_save_position || return;
        
        my %ignores;
        my ($field, $value);
        my $do_ignores = 0;
        while (my ($key, $val) = each %args) {
            if ($key =~ /^ignore_(\S+)/i) {
                my $ignore = uc($1);
                my $index  = $field_to_index{$ignore};
                $self->throw("'$ignore' wasn't a valid field name") unless defined $index;
                $ignores{$index} = $val;
                $do_ignores = 1;
            }
            else {
                ($field, $value) = ($key, $val);
            }
        }
        
        my $index;
        if ($field && $value) {
            $field = uc($field);
            $index = $field_to_index{$field};
        }
        
        my $pr = $self->parsed_record;
        my %lanes;
        
        # parse the whole file
        $self->_seek_first_record();
        RESULT: while ($self->next_record) {
            if ($do_ignores) {
                keys %ignores; # reset the iterator, since we may have nexted out
                while (my ($index, $regex) = each %ignores) {
                    next RESULT if $pr->[$index] =~ /$regex/i;
                }
            }
            
            if (defined $index) {
                if ($pr->[$index] eq $value) {
                    $lanes{ $pr->[2] } = 1;
                }
            }
            else {
                $lanes{ $pr->[2] } = 1;
            }
        }
        
        $self->_restore_position;
        
        return sort keys %lanes;
    }
    
    method _get_header {
        my $fh = $self->fh() || return;
        return 1 if $self->_header_parsed();
        
        my $saw = 0;
        while (<$fh>) {
            if (/^FASTQ_FILE\s+MD5\s+RUN_ID\s+STUDY_ID/) {
                $saw++;
                last;
            }
            else {
                # allow header line to not be present. lines can also end with
                # an extraneous \t\n, so remove that first... causes problems in
                # fake sequence.indexes when people leave the last few columns
                # empty...
                #s/\t\n$//;
                my @a = split("\t", $_);
                if (@a == 25 || @a == 26) {
                    my $tell = tell($fh);
                    if ($tell == -1) {
                        $self->throw("sequence.index has no header line, and you've piped it in - can't cope!");
                    }
                    seek($fh, 0, 0);
                    $saw++;
                    last;
                }
                else {
                    $self->warn("This file has no header line, and has " . scalar(@a) . " columns instead of 25 or 26");
                    last;
                }
            }
        }
        
        if ($saw) {
            $self->_set_header_parsed();
            return 1;
        }
        
        $self->throw("Unable to parse header before first result - is this a sequence.index file?");
    }
}

1;
