#!/usr/bin/env perl n 24 20:58:20 GMT 2013use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)]);
    use TestPipelines;
    use_ok('VRPipe::Steps::cufflinks');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('cufflinks', 'sam_files');
is_deeply [$step->id, $step->description], [1, 'Step for cufflinks: transcript assembly and quantification for rna-seq'], 'cufflinks step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'cufflinks',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        options => { delimiter => "\t" },
        source  => file(qw(t data cufflinks.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta => file(qw(t data sacCer3 chr1.fa))->absolute,
        #reference_fasta => file(qw(.. .. lust rna-seq general_files mouse_ref mm9.fa))->absolute,
        #known_genes_path => file(qw(.. .. lust rna-seq general_files knownGeneMm9.gtf))->absolute,
        #known_genes_path => file(qw(t data knownGeneMm9.gtf))->absolute,
        known_genes_path => file(qw(t data sacCer3 sacCer3_arabicChrA.chr1.gtf))->absolute,
        #gene_mask_path => file(qw(t data GeneMaskMm9.gtf))->absolute
    }
);

my @output_subdirs = output_subdirs(1);
my $outputfile_1   = file(@output_subdirs, '1_cufflinks', "transcripts.gtf");
my $outputfile_2   = file(@output_subdirs, '1_cufflinks', "isoforms.fpkm_tracking");
my $outputfile_3   = file(@output_subdirs, '1_cufflinks', "genes.fpkm_tracking");
my $outputfile_4   = file(@output_subdirs, '3_cufflinks', "skipped.gtf");
my @outputfiles;
push @outputfiles, $outputfile_1, $outputfile_2, $outputfile_3; # $outputfile_4; check these when I sort out better test data

#SKIP: {
#    skip "main test disabled due to lack of test data", 1;
ok handle_pipeline(@outputfiles), 'cufflinks pipeline ran ok, generating the expected files';
#}
