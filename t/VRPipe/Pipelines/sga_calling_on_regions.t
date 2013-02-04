#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(sga)]
    );
    use TestPipelines;
}

ok my $preprocess_pipeline = VRPipe::Pipeline->create(name => 'sga_prepare_region_fastq'), 'able to get the sga_prepare_region_fastq pipeline';
ok my $sga_pipeline        = VRPipe::Pipeline->create(name => 'sga_variant_calling'),      'able to get the sga_variant_calling pipeline';

my $calling_dir = get_output_dir('sga_calling_on_regions_test');

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz))->absolute);
my $ref_fa = file($calling_dir, 'human_g1k_v37.chr11.chr20.fa')->absolute->stringify;
my $new_ref_fa = VRPipe::File->create(path => $ref_fa);
my $oh         = $original_ref_fa->openr;
my $nh         = $new_ref_fa->openw;
while (<$oh>) {
    print $nh $_;
}
close($oh);
close($nh);

VRPipe::PipelineSetup->create(
    name       => 'sga prepare fastq test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data sga_calling_datasource.fofn))->absolute->stringify,
        options => { metadata_keys => 'sample|platform' },
    ),
    output_root => $calling_dir,
    pipeline    => $preprocess_pipeline,
    options     => {
        reference_index                   => file(qw(t data human_g1k_v37.chr11.chr20.fa.gz.fai))->absolute,
        split_bam_by_region_chunk_size    => 50000000,
        split_bam_by_region_chunk_overlap => 1000,
        split_bam_by_region_chrom_list    => '11 20',
        split_bam_by_region_include_mate  => 1,
        sga_preprocess_options            => '--min-length=50',
        sga_exe                           => 'sga',
        ignore_read_ordering              => 1,
        cleanup                           => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'sga calling test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '1[5]',
        options => { metadata_keys => 'continent|chrom|from|to' }
    ),
    output_root => $calling_dir,
    pipeline    => $sga_pipeline,
    options     => {
        reference_fasta   => $ref_fa,
        sga_exe           => 'sga',
        sga_index_options => '--no-reverse -t 4 -a ropebwt',
        cleanup           => 0
    }
);

ok handle_pipeline(), 'sga_prepare_fastq and sga_variant_calling pipelines ran ok';

my $chunks = [{ chrom => 11, from => 1, to => 50000000 }, { chrom => 11, from => 49999001, to => 99999000 }, { chrom => 11, from => 99998001, to => 135006516 }, { chrom => 20, from => 1, to => 50000000 }, { chrom => 20, from => 49999001, to => 63025520 }];

my @input_files;
my %samples = ('NA20340' => 1, 'HG02449' => 2, 'NA20281' => 3, 'HG01958' => 4, 'NA19381' => 5, 'NA19334' => 6);
while (my ($sample, $element_id) = each %samples) {
    my @output_subdirs = output_subdirs($element_id, 1);
    my $id = 0;
    foreach my $chunk (@$chunks) {
        my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
        push(@input_files, file(@output_subdirs, '1_bam_split_by_region', qq[$region.$sample.bam]));
        next unless ($$chunk{from} == 1); # not all chunks have reads
        push(@input_files, file(@output_subdirs, '3_bam_to_fastq',   qq[${region}.${sample}_$id.1.fastq]));
        push(@input_files, file(@output_subdirs, '3_bam_to_fastq',   qq[${region}.${sample}_$id.2.fastq]));
        push(@input_files, file(@output_subdirs, '4_sga_preprocess', qq[${region}.${sample}_$id.processed.fq.gz]));
        $id++;
    }
}

my @ref_files;
foreach my $suffix (qw(fa sai rsai bwt rbwt ssa)) {
    push @ref_files, file($calling_dir, 'human_g1k_v37.chr11.chr20.permute.' . $suffix);
}

my @calling_files;
foreach my $element_id (7, 8) {
    my @output_subdirs = output_subdirs($element_id, 2);
    foreach my $suffix (qw(fq.gz popidx bwt sai)) {
        push(@calling_files, file(@output_subdirs, '4_fastq_merge_and_index', 'merged.' . $suffix));
    }
    foreach my $suffix (qw(base.vcf calls.vcf evidence.bam variant.vcf)) {
        push(@calling_files, file(@output_subdirs, '6_sga_reference_based_calling', '0.merged.' . $suffix));
    }
}

ok handle_pipeline(@input_files, @ref_files, @calling_files), 'sga_prepare_fastq and sga_variant_calling pipelines created expected output files';

finish;
