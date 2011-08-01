use VRPipe::Base;

class VRPipe::Steps::bam_reheader with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                            max_files => -1,
                                                            description => '1 or more bam files',
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                         library => 'library name',
                                                                         sample => 'sample name',
                                                                         center_name => 'center name',
                                                                         platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                         study => 'name of the study',
                                                                         insert_size => 'expected (mean) insert size if paired',,
                                                                         bases => 'total number of base pairs',
                                                                         reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                         mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped',
                                                                         chunk => 'mapped_fastq(s) are this chunk of original fastq(s)',
                                                                         optional => ['library', 'insert_size', 'sample', 'center_name', 'platform', 'study']}),
                 dict_file => VRPipe::StepIODefinition->get(type => 'txt',
                                                            description => 'a sequence dictionary file for your reference fasta') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            $self->throw("bam_reheader is not yet implemented");
        };
    }
    method outputs_definition {
        return { headed_bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                          max_files => -1,
                                                          description => 'a bam file with good header',
                                                          metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                       library => 'library name',
                                                                       sample => 'sample name',
                                                                       center_name => 'center name',
                                                                       platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                       study => 'name of the study, put in the DS field of the RG header line',
                                                                       insert_size => 'expected (mean) insert size if paired',,
                                                                       bases => 'total number of base pairs',
                                                                       reads => 'total number of reads (sequences)',
                                                                       paired => '0=unpaired reads were mapped; 1=paired reads were mapped; 2=mixture of paired and unpaired reads were mapped',
                                                                       mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped',
                                                                       optional => ['library', 'insert_size', 'sample', 'center_name', 'platform', 'study']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Replaces a bam header so that it has complete sequence information, a good RG line, and chained PG lines";
    }
}

1;
