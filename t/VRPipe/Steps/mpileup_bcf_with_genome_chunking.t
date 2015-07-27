#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS HTSLIB)],
        required_exe => [qw(samtools bcftools htsfile)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::mpileup_bcf_with_genome_chunking');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('mpileup_bcf_with_genome_chunking', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Run mpileup for one or more BAM or CRAM files. Split calling into multiple jobs across the genome generating one BCF file per chunk.'], 'mpileup_bcf_with_genome_chunking step created and has correct description';

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz))->absolute);
my $ref_fa = file($output_dir, 'human_g1k_v37.chr11.chr20.fa')->absolute->stringify;
my $new_ref_fa = VRPipe::File->create(path => $ref_fa);
my $oh         = $original_ref_fa->openr;
my $nh         = $new_ref_fa->openw;
while (<$oh>) {
    print $nh $_;
}
close($oh);
close($nh);

system('samtools faidx ' . $ref_fa);

my $setup = VRPipe::PipelineSetup->create(
    name       => 'mpileup_bcf_with_genome_chunking',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'all',
        source  => file(qw(t data calling_datasource.fofn))->absolute->stringify,
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta          => $ref_fa,
        chrom_list               => '11 20',
        chunk_size               => 10000000,
        chunk_overlap            => 0,
        samtools_mpileup_options => '-EDSV -C50 -m2 -F0.0005 -d 10000 -g',
    }
);

ok handle_pipeline(), 'mpileup_bcf_with_genome_chunking pipeline ran ok';

finish;
