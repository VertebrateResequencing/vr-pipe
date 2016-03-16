#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(bcftools)]
    );
    use TestPipelines;
}

ok my $calling_pipeline = VRPipe::Pipeline->create(name => 'bam_calling_with_mpileup'), 'able to get the bam_calling_with_mpileup pipeline';
ok my $merge_pipeline   = VRPipe::Pipeline->create(name => 'bcftools_merge'),           'able to get the bcftools_merge pipeline';

my $calling_dir = get_output_dir('mpileup_gvcf');

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
    name       => 'mpileup gvcf test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'index' }
    ),
    output_root => $calling_dir,
    pipeline    => $calling_pipeline,
    options     => {
        samtools_exe             => '/nfs/users/nfs_p/pd3/git/pdtools/bcftools/bcftools-gvcf',
        samtools_mpileup_options => '-Ob -t AD,INFO/AD -C50 -pm2 -F0.1 -d10000 --gvcf 1,2,3,4,5,10,15',
        reference_fasta          => $ref_fa,
        chrom_list               => '11 20',
        chunk_size               => 0,
        chunk_overlap            => 0,
        bcftools_concat_options  => '-Ob -o $output_bcf',
        cleanup                  => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'mpileup merge gvcf test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_all',
        source  => '1[2]',
        options => {}
    ),
    output_root => $calling_dir,
    pipeline    => $merge_pipeline,
    options     => {
        bcftools_exe            => '/nfs/users/nfs_p/pd3/git/pdtools/bcftools/bcftools-gvcf',
        bcftools_merge_options  => '-Ou --gvcf | $bcftools call -Ob -mv -o $output_bcf',
        reference_fasta         => $ref_fa,
        chrom_list              => '11 20',
        chunk_size              => 0,
        chunk_overlap           => 0,
        bcftools_concat_options => '-Ob -o $output_bcf',
        cleanup                 => 0,
    }
);

ok handle_pipeline(), 'mpileup gvcf pipelines ran ok';

done_testing;
exit;
