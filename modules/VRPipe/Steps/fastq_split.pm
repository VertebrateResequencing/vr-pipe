use VRPipe::Base;

class VRPipe::Steps::fastq_split with VRPipe::StepRole {
    method options_definition {
        return { fastq_chunk_size => VRPipe::StepOption->get(description => 'when splitting fastq files into smaller chunks, this sets the size; defaults to something appropriate based on the mapper/platform', optional => 1) };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                              max_files => 3,
                                                              description => '1-3 fastq files',
                                                              metadata => {platform => 'the platform (aka technology), used to do the sequencing',
                                                                           reads => 'total number of reads (sequences)',
                                                                           paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                           mate => 'if paired, the path to the fastq that is our mate',
                                                                           optional => ['mate']}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
        };
    }
    method outputs_definition {
        return { split_fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                                    max_files => -1,
                                                                    description => 'split fastq files',
                                                                    metadata => {platform => 'the platform (aka technology), used to do the sequencing',
                                                                                 source_fastq => 'the fastq file this was split from',
                                                                                 reads => 'total number of reads (sequences) in this chunk',
                                                                                 paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                                 mate => 'if paired, the path to the fastq that is our mate and of the corresponding chunk',
                                                                                 optional => ['mate']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Takes a single-ended fastq file and/or 2 fastq files of a paired-end run and splits them into multiple smaller fastq files";
    }
}

1;
