#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES',
                    required_exe => [qw(bwa samtools)]);
    use TestPipelines;
}

my $mapping_output_dir = get_output_dir('bam_remapping');

ok my $mapping_pipeline = VRPipe::Pipeline->create(name => 'bam_mapping_with_bwa'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($mapping_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(sequence_dictionary bwa_index bam_metadata bam_name_sort bwa_aln_bam bwa_sam_using_bam sam_to_fixed_bam bam_reheader)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->create(name       => 's_suis remapping',
                                                          datasource => VRPipe::DataSource->create(type   => 'fofn',
                                                                                                   method => 'all',
                                                                                                   source => file(qw(t data datasource.bam_fofn))),
                                                          output_root => $mapping_output_dir,
                                                          pipeline    => $mapping_pipeline,
                                                          options     => {
                                                                       reference_fasta               => $ref_fa,
                                                                       reference_assembly_name       => 'SSuis1',
                                                                       reference_public_url          => 'ftp://s.suis.com/ref.fa',
                                                                       reference_species             => 'S.Suis',
                                                                       bwa_index_options             => '-a is',
                                                                       samtools_sort_options         => '-m 769M',
                                                                       uncompressed_fixed_bam_output => 0,
                                                                       cleanup                       => 0,
                                                                       sequence_dictionary_memory    => 150,
                                                                       sequence_dictionary_time      => 1,
                                                                       bwa_index_memory              => 150,
                                                                       bwa_index_time                => 1,
                                                                       bam_metadata_memory           => 150,
                                                                       bam_metadata_time             => 1,
                                                                       bwa_aln_bam_memory            => 150,
                                                                       bwa_aln_bam_time              => 1,
                                                                       bwa_sam_using_bam_memory      => 150,
                                                                       bwa_sam_using_bam_time        => 1,
                                                                       sam_to_fixed_bam_memory       => 150,
                                                                       sam_to_fixed_bam_time         => 1 });

ok handle_pipeline(), 'pipeline ran ok';

#*** needs proper tests...

finish;
