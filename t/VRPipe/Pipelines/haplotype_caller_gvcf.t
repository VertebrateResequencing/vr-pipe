#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK2 GATK_JAVA)],
        required_exe => [qw(bcftools tabix vcf-concat)]
    );
    use TestPipelines;
}

ok my $calling_pipeline        = VRPipe::Pipeline->create(name => 'bam_calling_with_gatk_haplotype_caller'),              'able to get the bam_calling_with_gatk_haplotype_caller pipeline';
ok my $gvcf_genotype_pipeline0 = VRPipe::Pipeline->create(name => 'gvcf_calling_with_gatk_haplotype_caller'),             'able to get the gatk_genotype_gvcfs pipeline';
ok my $gvcf_genotype_pipeline1 = VRPipe::Pipeline->create(name => 'gvcf_calling_one_combine_with_gatk_haplotype_caller'), 'able to get the gatk_genotype_gvcfs pipeline';

my $calling_dir = get_output_dir('haplotype_caller_gvcf_test');

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

copy(file(qw(t data human_g1k_v37.chr11.chr20.fa.gz.fai))->absolute->stringify,  file($calling_dir, 'human_g1k_v37.chr11.chr20.fa.fai')->absolute->stringify);
copy(file(qw(t data human_g1k_v37.chr11.chr20.fa.gz.dict))->absolute->stringify, file($calling_dir, 'human_g1k_v37.chr11.chr20.dict')->absolute->stringify);

my $override_file = file(qw(t data wgs_calling_override_options))->absolute->stringify;

my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data calling_datasource2.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'calling_datasource2.fofn'));
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

VRPipe::PipelineSetup->create(
    name       => 'haplotype caller gvcf test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'index' }
    ),
    output_root => $calling_dir,
    pipeline    => $calling_pipeline,
    options     => {
        gatk_path                => $ENV{GATK2},
        java_exe                 => $ENV{GATK_JAVA},
        reference_fasta          => $ref_fa,
        haplotype_caller_options => '--minPruning 3 --maxNumHaplotypesInPopulation 200 -ERC GVCF --max_alternate_alleles 3 -U ALLOW_SEQ_DICT_INCOMPATIBILITY -variant_index_type LINEAR -variant_index_parameter 128000',
        chrom_list               => '11 20',
        chunk_size               => 30000000,
        chunk_overlap            => 0,
        bcftools_concat_opts     => '',
        cleanup                  => 0,
    }
);

my $merge_setup1 = VRPipe::PipelineSetup->create(
    name       => 'GATK genotype gVCFs no combine',
    pipeline   => $gvcf_genotype_pipeline0,
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe_with_genome_chunking',
        method  => 'group_all',
        source  => '1[2]',
        options => {
            reference_index     => $ref_fa . '.fai',
            chunk_size          => 0,
            chrom_list          => '11 20',
            chunk_overlap       => 0,
            chunk_override_file => $override_file,
        },
    ),
    options => {
        reference_fasta => $ref_fa,
        chunk_size      => 50000000,
        chunk_overlap   => 0,
        chrom_list      => '11 20',
        gatk_path       => $ENV{GATK2},
        java_exe        => $ENV{GATK_JAVA},
        cleanup         => 0,
    },
    output_root => $calling_dir
);

my $merge_setup2 = VRPipe::PipelineSetup->create(
    name       => 'GATK genotype gVCFs one combine',
    pipeline   => $gvcf_genotype_pipeline1,
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe_with_genome_chunking',
        method  => 'group_all',
        source  => '1[2]',
        options => {
            reference_index     => $ref_fa . '.fai',
            chunk_size          => 0,
            chrom_list          => '11 20',
            chunk_override_file => $override_file,
        },
    ),
    options => {
        reference_fasta          => $ref_fa,
        maximum_gvcfs_to_combine => 2,
        chunk_size               => 50000000,
        chunk_overlap            => 0,
        chrom_list               => '11 20',
        gatk_path                => $ENV{GATK2},
        java_exe                 => $ENV{GATK_JAVA},
        cleanup                  => 0,
    },
    output_root => $calling_dir
);

ok handle_pipeline(), 'gvcf pipelines ran ok';

done_testing;
exit;
