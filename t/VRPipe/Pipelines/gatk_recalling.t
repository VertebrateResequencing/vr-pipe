#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK)],
        required_exe => [qw(samtools bcftools tabix)],
    );
    use TestPipelines;
}

ok my $mpileup_pipeline = VRPipe::Pipeline->get(name => 'snp_calling_mpileup'), 'able to get the snp_calling_mpileup pipeline';
ok my $gatk_pipeline = VRPipe::Pipeline->get(name => 'snp_calling_gatk_unified_genotyper'), 'able to get the snp_calling_gatk_unified_genotyper pipeline';

my $calling_dir = get_output_dir('gatk_recalling_test');

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

my $all_samples_file = file(qw(t data wgs_calling.samples))->absolute->stringify;

my $original_override_file = VRPipe::File->create(path => file(qw(t data wgs_calling_override_options))->absolute);
my $override_file = VRPipe::File->create(path => file($calling_dir, 'wgs_calling_override_options')->absolute);
$override_file->touch;

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data wgs_calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'wgs_calling_datasource.fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
print $ofh scalar <$ifh>;
while (<$ifh>) {
    chomp;
    my ($source_path, @meta) = split(/\t/, $_);
    my $source = file($source_path);
    my $dest = file($calling_dir, $source->basename);
    copy($source, $dest);
    print $ofh join("\t", $dest, @meta);
    print $ofh "\n";
}
$orig_fofn_file->close;
$fofn_file->close;

VRPipe::PipelineSetup->create(
    name       => 'wgs test multi-sample mpileup calling',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata_with_genome_chunking',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path,
        options => {
            metadata_keys       => 'analysis_group|chrom',
            reference_index     => file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify,
            chrom_list          => '11 20',
            chunk_size          => 10000000,
            chunk_overlap       => 0,
            chunk_override_file => $override_file->path->stringify
        },
    ),
    output_root => $calling_dir,
    pipeline    => $mpileup_pipeline,
    options     => {
        reference_fasta          => $ref_fa,
        samtools_mpileup_options => '-EDSV -C50 -m2 -F0.0005 -d 10000 -ug',
        cleanup                  => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'gatk recalling',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => '1[0,2]', # BAMs and VCFs
        options => {},
    ),
    output_root => $calling_dir,
    pipeline    => $gatk_pipeline,
    options     => {
        reference_assembly_name   => 'NCBI37',
        reference_species         => 'Human',
        reference_fasta           => $ref_fa,
        unified_genotyper_options => '--genotype_likelihoods_model BOTH -stand_call_conf 0.0 -stand_emit_conf 0.0 -baq RECALCULATE -out_mode EMIT_ALL_SITES',
        cleanup                   => 0,
    }
);

my $chunks = [{ chrom => 11, from => 1, to => 10000000 }, { chrom => 11, from => 10000001, to => 20000000 }, { chrom => 11, from => 20000001, to => 30000000 }, { chrom => 11, from => 30000001, to => 40000000 }, { chrom => 11, from => 40000001, to => 50000000 }, { chrom => 11, from => 50000001, to => 60000000 }, { chrom => 11, from => 60000001, to => 70000000 }, { chrom => 11, from => 70000001, to => 80000000 }, { chrom => 11, from => 80000001, to => 90000000 }, { chrom => 11, from => 90000001, to => 100000000 }, { chrom => 11, from => 100000001, to => 110000000 }, { chrom => 11, from => 110000001, to => 120000000 }, { chrom => 11, from => 120000001, to => 130000000 }, { chrom => 11, from => 130000001, to => 135006516 }, { chrom => 20, from => 1, to => 10000000 }, { chrom => 20, from => 10000001, to => 20000000 }, { chrom => 20, from => 20000001, to => 30000000 }, { chrom => 20, from => 30000001, to => 40000000 }, { chrom => 20, from => 40000001, to => 50000000 }, { chrom => 20, from => 50000001, to => 60000000 }, { chrom => 20, from => 60000001, to => 63025520 }];

my $element_id = 1;
my @calling_files;
foreach my $chunk (@$chunks) {
    my @output_subdirs = output_subdirs($element_id, 1);
    my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
    push(@calling_files, file(@output_subdirs, '2_mpileup_vcf', qq[$region.mpileup.vcf.gz]));
    push(@calling_files, file(@output_subdirs, '2_mpileup_vcf', qq[$region.mpileup.vcf.gz.tbi]));
    ++$element_id;
}

ok handle_pipeline(@calling_files), 'pipeline ran ok and all calling files were created';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->get(pipelinesetup => 2, stepmember => 7, dataelement => $element_id)->cmd_summary->summary], ['samtools mpileup -EDSV -C50 -m2 -F0.0005 -d 10000 -ug -r $region -f $reference_fasta -b $bams_list | bcftools view -p 0.99 -vcgN -s $samples_file - | bgzip -c > $vcf_file', 'java $jvm_args -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference_fasta -I $bams_list -o $vcf_file --genotype_likelihoods_model BOTH -stand_call_conf 0.0 -stand_emit_conf 0.0 -baq RECALCULATE -out_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles $sites_file -L $region'], 'cmd summaries for the major steps were as expected';

finish;
