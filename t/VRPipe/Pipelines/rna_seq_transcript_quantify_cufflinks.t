#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES CUFFLINKS_EXE CUFFLINKS_GENOME_FASTA_PATH CUFFLINKS_KNOWN_GENES_PATH CUFFLINKS_GENE_MASK_PATH)],
        required_exe => []
    );
    use TestPipelines;
}
my $output_dir = get_output_dir('rna_seq_transcript_quantify_cufflinks-test');
ok my $pipeline = VRPipe::Pipeline->create(name => 'rna_seq_transcript_quantify_cufflinks'), 'able to get the rna_seq_transcript_quantify_cufflinks pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}

is_deeply \@s_names, [qw(cufflinks)], 'the pipeline has the correct steps';

my $pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'rna_seq_transcript_quantify_cufflinks',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        options => { delimiter => "\t" },
        source  => file(qw(t data cufflinks.fofn)),
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
);
my @output_subdirs = output_subdirs(1);

my $outputfile_1 = file(@output_subdirs, '1_cufflinks', "transcripts.gtf");
my @outputfiles;
push(@outputfiles, $outputfile_1);
ok handle_pipeline(@outputfiles), 'rna-seq-pipeline ran ok, generating the expected output files';
finish;
