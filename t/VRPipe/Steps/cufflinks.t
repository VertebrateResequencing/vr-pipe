#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES CUFFLINKS_EXE CUFFLINKS_GENOME_FASTA_PATH CUFFLINKS_KNOWN_GENES_PATH CUFFLINKS_GENE_MASK_PATH)]);
    use TestPipelines;
    use_ok('VRPipe::Steps::cufflinks');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('cufflinks', 'sam_files');
is_deeply [$step->id, $step->description], [1, 'Step for cufflinks: transcript assembly and quantification for rna-seq'], 'cufflinks step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(name       => 'cufflinks',
                                          datasource => VRPipe::DataSource->create(type    => 'delimited',
                                                                                   method  => 'all_columns',
                                                                                   options => { delimiter => "\t" },
                                                                                   source  => file(qw(t data cufflinks.fofn))->absolute),
                                          output_root => $output_dir,
                                          pipeline    => $pipeline,
                                          options     => {});

my @output_subdirs = output_subdirs(1);
my $outputfile_1   = file(@output_subdirs, '1_cufflinks', "transcripts.gtf");
my $outputfile_2   = file(@output_subdirs, '1_cufflinks', "isoforms.fpkm_tracking");
my $outputfile_3   = file(@output_subdirs, '1_cufflinks', "genes.fpkm_tracking");
my @outputfiles;
push @outputfiles, $outputfile_1, $outputfile_2, $outputfile_3; # $outputfile_4;
ok handle_pipeline(@outputfiles), 'cufflinks pipeline ran ok, generating the expected file';

