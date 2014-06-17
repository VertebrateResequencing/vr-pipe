#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(bamcheck bam2fastq samtools STAR)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('rna_seq_star_bam_remapping');

# check pipeline has correct steps
ok my $star_pipeline = VRPipe::Pipeline->create(name => 'rna_seq_star_bam_remapping'), 'able to get the rna_seq_star_bam_remapping pipeline';
my @sb_names;
foreach my $stepmember ($star_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(sequence_dictionary star_buildgenome bam_metadata bam_to_fastq star_map_fastq sam_to_fixed_bam bam_reheader bam_index)], 'the rna_seq_star_bam_remapping pipeline has the correct steps';

my $ref_fasta_orig = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($output_dir, 'ref');
$star_pipeline->make_path($ref_dir);
my $ref_fasta = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fasta_orig, $ref_fasta);

my $star_setup = VRPipe::PipelineSetup->create(
    name       => 'star setup',
    pipeline   => $star_pipeline,
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data hs_chr20.qc.bam.fofn)),
    ),
    options => {
        reference_fasta => $ref_fasta,
    },
    output_root => $output_dir
);

my @final_files;
foreach my $name (qw(Genome SA SAindex)) {
    push(@final_files, file($ref_dir, $name));
}
foreach my $element (@{ get_elements($star_setup->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $star_setup->id);
    push(@final_files, file(@output_dirs, '7_bam_reheader', 'Aligned.out.bam.bai'));
}
ok handle_pipeline(@final_files), 'rna_seq_star_bam_remapping pipeline ran ok';

finish;
exit;
