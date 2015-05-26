#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK2)],
        required_exe => [qw(samtools bcftools htsfile)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::gatk_haplotype_caller_with_genome_chunking');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gatk_haplotype_caller_with_genome_chunking', 'bam_files', 'bai_files');
is_deeply [$step->id, $step->description], [1, 'Run GATK HaplotypeCaller for one or more BAMs, generating one compressed VCF per set of BAMs. Sites list can be provided'], 'gatk_haplotype_caller_with_genome_chunking step created and has correct description';

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

my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($output_dir, 'calling_datasource.fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
my $idx = 0;
while (<$ifh>) {
    chomp;
    my @a = split /\t/;
    $a[0] = $a[0] . ".bai";
    if ($idx) {
        print $ofh "$_\t$idx\n";
        print $ofh join("\t", @a, "$idx\n");
    }
    else {
        print $ofh "$_\tindex\n";
    }
    $idx++;
}
$orig_fofn_file->close;
$fofn_file->close;

system("samtools faidx $ref_fa");
my $dict = file($output_dir, 'human_g1k_v37.chr11.chr20.dict')->absolute->stringify;
system("samtools dict $ref_fa > $dict");

my $setup = VRPipe::PipelineSetup->create(
    name       => 'gatk_haplotype_caller_with_genome_chunking',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'index' }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        reference_fasta          => $ref_fa,
        chrom_list               => '11 20',
        chunk_size               => 10000000,
        chunk_overlap            => 0,
        haplotype_caller_options => '--genotyping_mode DISCOVERY -stand_call_conf 0.0 -stand_emit_conf 0.0',
    }
);

ok handle_pipeline(), 'gatk_haplotype_caller_with_genome_chunking pipeline ran ok';

finish;
