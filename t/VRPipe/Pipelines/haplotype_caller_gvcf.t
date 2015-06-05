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

ok my $calling_pipeline       = VRPipe::Pipeline->create(name => 'sample_gvcf_with_gatk_haplotype_caller'), 'able to get the gvcf_gatk_haplotype_caller pipeline';
ok my $gvcf_genotype_pipeline = VRPipe::Pipeline->create(name => 'gatk_genotype_gvcfs'),                    'able to get the gatk_genotype_gvcfs pipeline';
ok my $concat_pipeline        = VRPipe::Pipeline->create(name => 'vcf_concat'),                             'able to get the vcf_concat pipeline';

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

my $override_file = file(qw(t data wgs_calling_override_options))->absolute->stringify;

my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'calling_datasource.fofn'));
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
        reference_assembly_name  => 'NCBI37',
        reference_species        => 'Human',
        reference_fasta          => $ref_fa,
        haplotype_caller_options => '--minPruning 3 --maxNumHaplotypesInPopulation 200 -ERC GVCF --max_alternate_alleles 3 -U ALLOW_SEQ_DICT_INCOMPATIBILITY -variant_index_type LINEAR -variant_index_parameter 128000',
        chrom_list               => '11 20',
        chunk_size               => 10000000,
        chunk_overlap            => 0,
        bcftools_concat_opts     => '',
        cleanup                  => 0,
    }
);

my $merge_setup = VRPipe::PipelineSetup->create(
    name       => 'GATK genotype gVCFs',
    pipeline   => $gvcf_genotype_pipeline,
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe_with_genome_chunking',
        method  => 'group_all',
        source  => '1[4]',
        options => {
            reference_index     => $ref_fa . '.fai',
            chunk_size          => 10000000,
            chunk_overlap       => 0,
            chunk_override_file => $override_file,
        },
    ),
    options => {
        reference_fasta          => $ref_fa,
        maximum_gvcfs_to_combine => 3,
        gatk_path                => $ENV{GATK2},
        java_exe                 => $ENV{GATK_JAVA},
    },
    output_root => $calling_dir
);

VRPipe::PipelineSetup->create(
    name       => 'vcf-concat test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '2[2]',
        options => { metadata_keys => 'chrom' }
    ),
    output_root => $calling_dir,
    pipeline    => $concat_pipeline,
    options     => {}
);

ok handle_pipeline(), 'gvcf pipelines ran ok';

done_testing;
exit;
