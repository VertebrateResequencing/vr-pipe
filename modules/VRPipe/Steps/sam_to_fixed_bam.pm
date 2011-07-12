use VRPipe::Base;

class VRPipe::Steps::sam_to_fixed_bam with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { sam_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                            max_files => 2,
                                                            description => '1-2 sam files',
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run)',
                                                                         reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                         mate => 'if paired, the path to the fastq that is our mate',
                                                                         optional => ['mate']}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("foo");
        };
    }
    method outputs_definition {
        return { fixed_bam_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                                  max_files => 2,
                                                                  description => 'mapped bam file(s)',
                                                                  metadata => {lane => 'lane name (a unique identifer for this sequencing run)',
                                                                               reads => 'total number of reads (sequences)',
                                                                               paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                               mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Turns a sam file into a well-constructed bam file, with full header, sorting, fixed mates, correct NM values and RG tagged on every record";
    }
}

1;
