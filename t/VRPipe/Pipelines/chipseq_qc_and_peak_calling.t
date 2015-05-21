#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_PHANTOMPEAK_RSCRIPT VRPIPE_PHANTOMPEAK_RLIBS)],
        required_exe => [qw(samtools bedtools bamcheck macs2 bedGraphToBigWig Rscript)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('chipseq_qc_and_peak_calling');

# check pipeline has correct steps
ok my $chipseq_pipeline = VRPipe::Pipeline->create(name => 'chipseq_qc_and_peak_calling'), 'able to get the chipseq_qc_and_peak_calling pipeline';
my @sb_names;
foreach my $stepmember ($chipseq_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(bam_processing bam_metadata phantompeakqualtools macs_callpeak bedgraph2bigwig)], 'the chipseq_qc_and_peak_calling pipeline has the correct steps';

my $chipseq_setup = VRPipe::PipelineSetup->create(
    name       => 'chipseq setup',
    pipeline   => $chipseq_pipeline,
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        options => { metadata_keys => 'individual' },
        source  => file(qw(t data chipseq-bams.fofn)),
    ),
    options => {
        r_bin_path             => "Rscript",
        r_libs                 => $ENV{VRPIPE_PHANTOMPEAK_RLIBS},
        phantompeak_script     => $ENV{VRPIPE_PHANTOMPEAK_RSCRIPT},
        input_metadata_to_keep => 'individual,sample',
        sample_metadata_key    => 'sample',
        command_line           => q[samtools-exp-rc view -h -u -f 3 -F 12 -q 10 $input_bam | samtools-exp-rc sort -n -Osam -T samtools_nsort_tmp - | awk '/^@/; /\tX0:i:1\t/ { if ($1==qn){print l"\n"$0}; l=$0; qn=$1;}' | samtools-exp-rc sort -Obam -T samtools_csort_tmp - > $output_bam],
        reference_index        => file(qw(t data hs37d5.fa.fai))->absolute,
        macs2_callpeak_options => "K27AC:-Bg hs --fix-bimodal,K4ME3:-Bg hs --fix-bimodal,K27ME3:-Bg hs --broad --fix-bimodal,INPUT:-Bg hs --broad --fix-bimodal"
    },
    output_root => $output_dir
);

my @final_files;
my @sample_names = qw(coxy33_input coxy33_k27ac coxy33_k27me3 coxy33_k4me3);
foreach my $sample (@sample_names) {
    push(@final_files, file(output_subdirs(1), '5_bedgraph2bigwig', "${sample}_control_lambda.bdg.bw"));
    push(@final_files, file(output_subdirs(1), '5_bedgraph2bigwig', "${sample}_treat_pileup.bdg.bw"));
}
ok handle_pipeline(@final_files), 'chipseq_qc_and_peak_calling pipeline ran ok';

finish;
exit;
