use VRPipe::Base;

class VRPipe::Parser::bamcheck with VRPipe::ParserRole {
    has 'sequences' => (is => 'ro',
                        isa => 'Int',
                        writer => '_sequences');
    
    has 'is_paired' => (is => 'ro',
                        isa => 'Bool',
                        writer => '_is_paired');
    
    has 'is_sorted' => (is => 'ro',
                        isa => 'Bool',
                        writer => '_is_sorted');
    
    has 'first_fragments' => (is => 'ro',
                              isa => 'Int',
                              writer => '_first_fragments');
    
    has 'last_fragments' => (is => 'ro',
                             isa => 'Int',
                             writer => '_last_fragments');
    
    has 'reads_mapped' => (is => 'ro',
                           isa => 'Int',
                           writer => '_reads_mapped');
    
    has 'reads_unmapped' => (is => 'ro',
                           isa => 'Int',
                           writer => '_reads_unmapped');
    
    has 'reads_unpaired' => (is => 'ro',
                           isa => 'Int',
                           writer => '_reads_unpaired');
    
    has 'reads_paired' => (is => 'ro',
                           isa => 'Int',
                           writer => '_reads_paired');
    
    has 'reads_duplicated' => (is => 'ro',
                           isa => 'Int',
                           writer => '_reads_duplicated');
    
    has 'reads_mq0' => (is => 'ro',
                           isa => 'Int',
                           writer => '_reads_mq0');
    
    has 'total_length' => (is => 'ro',
                           isa => 'Int',
                           writer => '_total_length');
    
    has 'bases_mapped' => (is => 'ro',
                           isa => 'Int',
                           writer => '_bases_mapped');
    
    has 'bases_mapped_cigar' => (is => 'ro',
                           isa => 'Int',
                           writer => '_bases_mapped_cigar');
    
    has 'bases_trimmed' => (is => 'ro',
                           isa => 'Int',
                           writer => '_bases_trimmed');
    
    has 'bases_duplicated' => (is => 'ro',
                           isa => 'Int',
                           writer => '_bases_duplicated');
    
    has 'mismatches' => (is => 'ro',
                           isa => 'Int',
                           writer => '_mismatches');
    
    has 'error_rate' => (is => 'ro',
                           isa => 'Num',
                           writer => '_error_rate');
    
    has 'average_length' => (is => 'ro',
                         isa => 'Num',
                         writer => '_average_length');
    
    has 'maximum_length' => (is => 'ro',
                         isa => 'Int',
                         writer => '_maximum_length');
    
    has 'average_quality' => (is => 'ro',
                         isa => 'Num',
                         writer => '_average_quality');
    
    has 'insert_size_average' => (is => 'ro',
                         isa => 'Num',
                         writer => '_insert_size_average');
    
    has 'insert_size_standard_deviation' => (is => 'ro',
                                  isa => 'Num',
                                  writer => '_insert_size_standard_deviation');
    
    has 'inward_oriented_pairs' => (is => 'ro',
                         isa => 'Int',
                         writer => '_inward_oriented_pairs');
    
    has 'outward_oriented_pairs' => (is => 'ro',
                         isa => 'Int',
                         writer => '_outward_oriented_pairs');
    
    has 'pairs_with_other_orientation' => (is => 'ro',
                         isa => 'Int',
                         writer => '_pairs_with_other_orientation');
    
=head2 parsed_record

 Title   : parsed_record
 Usage   : my $parsed_record= $obj->parsed_record()
 Function: Get the data structure that will hold the last parsed record
           requested by next_record()
 Returns : array ref, currently empty, since only SN lines are parsed so far,
           and their content can be found with other methods
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
        
        # *** actual parsing of the data not yet supported
        return 0;
        
        return 1;
    }
    
    method _get_header {
        my $fh = $self->fh() || return;
        return 1 if $self->_header_parsed();
        
        #SN      sequences:      2000
        #SN      is paired:      1
        #SN      is sorted:      1
        #SN      1st fragments:  1000
        #SN      last fragments: 1000
        #SN      reads mapped:   1084
        #SN      reads unmapped: 916
        #SN      reads unpaired: 14
        #SN      reads paired:   1070
        #SN      reads duplicated:       0
        #SN      total length:   115000
        #SN      bases mapped:   62288
        #SN      bases mapped (cigar):   58663
        #SN      bases trimmed:  0
        #SN      bases duplicated:       0
        #SN      mismatches:     1262
        #SN      error rate:     2.151271e-02
        #SN      average length: 57
        #SN      maximum length: 61
        #SN      average quality:        20.9
        #SN      insert size average:    284.1
        #SN      insert size standard deviation: 70.9
        
        my $saw = 0;
        while (<$fh>) {
            if (/^SN\s+([^:]+):\s+(\S+)/) {
                my $method = $1;
                my $value = $2;
                my $orig_method = $method;
                $method =~ s/\s+/_/g;
                $method = 'first_fragments' if $method eq '1st_fragments';
                $method = 'bases_mapped_cigar' if $method eq 'bases_mapped_(cigar)';
                $method = lc('_'.$method);
                unless ($self->can($method)) {
                    $self->warn("unexpected SN line $orig_method");
                    next;
                }
                $self->$method($value);
                $saw++;
            }
            else {
                last if $saw >= 22;
            }
        }
        
        if ($saw >= 22) {
            $self->_set_header_parsed();
            return 1;
        }
        
        $self->throw("Unable to parse all SN lines (only saw $saw) - is this a bamcheck file?")
    }
}

1;