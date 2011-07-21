use VRPipe::Base;

class VRPipe::Steps::sam_to_fixed_bam with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { sam_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                            max_files => -1,
                                                            description => 'raw sam files from a mapper',
                                                            metadata => {reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired reads were mapped; 1=paired reads were mapped'}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("foo");
        };
    }
    method outputs_definition {
        return { fixed_bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                                  max_files => -1,
                                                                  description => 'uncompressed coordinate-sorted bam file(s)',
                                                                  metadata => {reads => 'total number of reads (sequences)',
                                                                               paired => '0=unpaired reads were mapped; 1=paired reads were mapped'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Turns a sam file into an uncompressed coordinate-sorted bam file with fixed mates and correct NM tag values";
    }
}

1;

=pod
 samtools view -bS
 samtools sort -n
 samtools fixmate
 samtools sort
 samtools fillmd -b
=cut