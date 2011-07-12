use VRPipe::Base;

class VRPipe::Steps::bwa_sam with VRPipe::StepRole {
    method options_definition {
        return { bwa_samse_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa samse command line, including desired options, but excluding the input sai and fastq files',
                                                          optional => 1,
                                                          default_value => 'bwa samse'),
                 bwa_sampe_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa sampe command line, including desired options, but excluding the input sai and fastq files; defaults to bwa sampe with -a set as appropriate for the fastq insert_size',
                                                          optional => 1) };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                              max_files => 3,
                                                              description => '1-3 fastq files',
                                                              metadata => {lane => 'lane name (a unique identifer for this sequencing run)',
                                                                           insert_size => 'expected (mean) insert size; specify as 0 if all fastqs are single ended',
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
        return { bwa_sam_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                                max_files => 2,
                                                                description => 'mapped sam file(s)',
                                                                metadata => {lane => 'lane name (a unique identifer for this sequencing run)',
                                                                             reads => 'total number of reads (sequences)',
                                                                             paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                             mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Produces sam files with bwa samse/sampe";
    }
}

1;
