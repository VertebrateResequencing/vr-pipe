#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(bamtofastq kallisto samtools)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('rna_seq_kallisto_quantification');

# check pipeline has correct steps
ok my $kallisto_pipeline = VRPipe::Pipeline->create(name => 'rna_seq_kallisto_quantification'), 'able to get the rna_seq_kallisto_quantification pipeline';
my @sb_names;
foreach my $stepmember ($kallisto_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(bamtofastq kallisto_index kallisto_quant)], 'the rna_seq_kallisto_quantification pipeline has the correct steps';

my $trans_fasta_orig = file(qw(t data transcripts.fasta.gz));
my $trans_dir = dir($output_dir, 'trans');
$kallisto_pipeline->make_path($trans_dir);
my $transcripts_fasta = file($trans_dir, 'transcripts.fasta.gz')->stringify;
copy($trans_fasta_orig, $transcripts_fasta);

my $kallisto_setup = VRPipe::PipelineSetup->create(
    name       => 'kallisto setup',
    pipeline   => $kallisto_pipeline,
    datasource => VRPipe::DataSource->create(
        type   => 'fofn_with_metadata',
        method => 'all',
        source => file(qw(t data reads.rnaseq.bam.fofn)),
    ),
    options => {
        kallisto_h5dump        => 1,
        transcripts_fasta      => $transcripts_fasta,
        kallisto_quant_options => "-b 10 --pseudobam",
    },
    output_root => $output_dir
);

my @final_files;
push(@final_files, file($trans_dir, "transcripts.fasta.gz_kallisto_index"));
foreach my $element (@{ get_elements($kallisto_setup->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $kallisto_setup->id);
    push(@final_files, file(@output_dirs, '3_kallisto_quant', 'pseudo_align.bam'));
    for my $i (0 .. 9) {
        push(@final_files, file(@output_dirs, '3_kallisto_quant', 'h5dump', "bs_abundance_$i.tsv"));
    }
}
ok handle_pipeline(@final_files), 'rna_seq_kallisto_bam_remapping pipeline ran ok';

finish;
exit;
