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

=pod
@SQ     SN:GL000192.1   LN:547496       UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  AS:NCBI37       M5:325ba9e808f669dfeee210fdd7b470ac     SP:Human
@RG     ID:1    PL:ILLUMINA     PU:110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9        LB:HUMbpgRUUDIAAPEI-9   PI:476  DS:Study UK10K_COHORT_TWINSUK   DT:2010-12-15T00:00:00+0000     SM:UK10K_QTL211459      CN:SC
@PG     ID:filter_adapter       PN:filter_adapter.pl    VN:1.0  DS:filter adapters      CL:perl /software/filter_adapter.pl  /HKC10040_HUMbpgR/HUMbpgRUUDIAAPEI-9/110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9/read_1.adapter.list /HKC10040_HUMbpgR/HUMbpgRUUDIAAPEI-9/110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9/read_2.adapter.list /HKC10040_HUMbpgR/HUMbpgRUUDIAAPEI-9/110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9/110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9_1.fq.gz /HKC10040_HUMbpgR/HUMbpgRUUDIAAPEI-9/110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9/110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9_2.fq.gz read_1.fq.gz read_2.fq.gz
@PG     ID:bwa_aln      PN:bwa  PP:filter_adapter       VN:0.5.9rc1 (r1561)     DS:bwa alignment from fastq files       CL:/software/bwa aln -q 15 -t 4 /data/human_g1k_v37.fasta read_1.fq.gz > 1.sai;/software/bwa aln -q 15 -t 4 /data/human_g1k_v37.fasta read_2.fq.gz > 2.sai
@PG     ID:bwa_sam      PN:bwa  PP:bwa_aln      VN:0.5.9rc1 (r1561)     DS:bwa converting alignemtns to sam     CL:/software/bwa sampe /data/human_g1k_v37.fasta 1.sai 2.sai read_1.fq.gz read_2.fq.gz > 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sam
@PG     ID:sam_to_bam   PN:samtools     PP:bwa_sam      VN:0.1.8 (r613) DS:convert sam to bam   CL:/software/samtools view -b -u -S 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sam > 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.bam
@PG     ID:bam_sort     PN:samtools     PP:sam_to_bam   VN:0.1.8 (r613) DS:sort bam file        CL:/software/samtools sort 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.bam 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sort
@PG     ID:remove_duplication   PN:samtools     PP:bam_sort     VN:0.1.8 (r613) DS:remove duplication   CL:/software/samtools rmdup 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sort.bam 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sort.rmdup.bam
@PG     ID:add_header   PN:samtools     PP:remove_duplication   VN:0.1.8 (r613) DS:add PG, RG and SQ tag in bam header  CL:/software/samtools reheader bam.header 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sort.rmdup.bam > 110310_I635_FC81DAEABXX_L7_HUMbpgRUUDIAAPEI-9.sort.rmdup.reheader.bam


@PG     ID:bwa_aln      PN:bwa aln      PP:fastq_creator        VN:0.5.8c (r1536)       DS:bwa alignment from fastq files       CL:/software/solexa/bin/aligners/bwa/bwa-0.5.8c/bwa aln -q 15 -t 2 /nfs/repository/d0031/references/Homo_sapiens/1000Genomes/all/bwa/human_g1k_v37.fasta /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_1.fastq > /tmp/sF99t1_QNN/1.sai;/software/solexa/bin/aligners/bwa/bwa-0.5.8c/bwa aln -q 15 -t 2 /nfs/repository/d0031/references/Homo_sapiens/1000Genomes/all/bwa/human_g1k_v37.fasta /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_2.fastq > /tmp/sF99t1_QNN/2.sai;
@PG     ID:bwa_sam      PN:bwa sam      PP:bwa_aln      VN:0.5.8c (r1536)       DS:bwa converting alignemtns to sam     CL:/software/solexa/bin/aligners/bwa/bwa-0.5.8c/bwa sampe /nfs/repository/d0031/references/Homo_sapiens/1000Genomes/all/bwa/human_g1k_v37.fasta /tmp/sF99t1_QNN/1.sai /tmp/sF99t1_QNN/2.sai /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_1.fastq /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_2.fastq
@PG     ID:sam_fastq_check      PN:sam_fastq_check      PP:bwa_sam      DS:program used to check bases and quality scores in bam are the same as the ones in original fastq files used to generate the bam      CL:/software/solexa/bin/sam_fastq_check.pl --fastq_files /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_1.fastq --fastq_files /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_2.fastq
@PG     ID:sam_second_call      PN:sam_second_call      PP:sam_fastq_check      VN:11838        DS:add second call bases to each read as E2 string      CL:/software/solexa/bin/sam_second_call.pl --lane 8 --qseq_dir /nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/SecondCall --read_number_1 1 --read_number_2 2
@PG     ID:sam_header   PN:sam_header   PP:sam_second_call      VN:12169        DS:add PG, RG and SQ tag in bam header, and add RG tag for each read
@PG     ID:picard_samformatconverter    PN:Picard SamFormatConverter    PP:sam_header   VN:1.39 DS:Picard used to convert sam to bam    CL:/software/bin/java -jar /software/solexa/bin/aligners/picard/picard-tools-1.39/SamFormatConverter.jar VALIDATION_STRINGENCY=SILENT INPUT=/dev/stdin OUTPUT=/nfs/sf39/ILorHSany_sf39/analysis/110206_IL22_05794/Data/Intensities/Bustard1.8.1a2_14-02-2011_RTA/GERALD_14-02-2011_RTA/archive/5794_8_unsorted.bam
@PG     ID:samtools_sort        PN:samtools     PP:picard_samformatconverter    VN:0.1.12a (r862)       DS:Samtools used to sort reads in bam   CL:/software/solexa/bin/aligners/samtools/samtools-0.1.12a/samtools sort
@PG     ID:picard_markduplicates        PN:Picard MarkDuplicates        PP:samtools_sort        VN:1.39 DS:Picard used to mark duplicates reads in bam file     CL:/software/bin/java -Xmx3000m -jar /software/solexa/bin/aligners/picard/picard-tools-1.39/MarkDuplicates.jar INPUT= OUTPUT= METRICS_FILE= TMP_DIR= VALIDATION_STRINGENCY='SILENT' MAX_FILE_HANDLES_FOR_READ_ENDS_MAP='900' REMOVE_DUPLICATES='false' ASSUME_SORTED='true' VERBOSITY='ERROR' READ_NAME_REGEX='[a-zA-Z0-9_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*' CREATE_MD5_FILE='TRUE'
=cut
