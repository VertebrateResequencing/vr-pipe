#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK)]);
    use TestPipelines;
}

ok my $calling_pipeline = VRPipe::Pipeline->get(name => 'snp_calling_gatk_unified_genotyper'), 'able to get the snp_calling_gatk_unified_genotyper pipeline';
ok my $concat_pipeline = VRPipe::Pipeline->get(name => 'vcf_concat'), 'able to get the vcf_concat pipeline';

my $calling_dir = get_output_dir('wgs_single_sample_gatk_calling_test');

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data S_suis_P17.fa))->absolute);
my $ref_fa = VRPipe::File->create(path => file($calling_dir, 'S_suis_P17.fa')->absolute);
$original_ref_fa->copy($ref_fa);

my $override_file = file(qw(t data wgs_calling_override_options))->absolute->stringify;

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data wgs_single_sample_calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'wgs_single_sample_calling_datasource.fofn'));
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
    name       => 'wgs test single-sample calling',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata_with_genome_chunking',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path,
        options => {
            reference_index     => file(qw(t data S_suis_P17.fa.fai))->absolute->stringify,
            chrom_list          => 'fake_chr1 fake_chr2',
            chunk_size          => 100000,
            chunk_overlap       => 0,
            chunk_override_file => $override_file,
            metadata_keys       => 'sample'
        }
    ),
    output_root => $calling_dir,
    pipeline    => $calling_pipeline,
    options     => {
        reference_fasta           => $ref_fa->path,
        unified_genotyper_options => '--genotype_likelihoods_model BOTH -stand_call_conf 0.0 -stand_emit_conf 0.0',
        cleanup                   => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'vcf-concat test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '1[4]',
        options => { metadata_keys => 'sample' }
    ),
    output_root => $calling_dir,
    pipeline    => $concat_pipeline,
    options     => {}
);

my $chunks = [{ chrom => 'fake_chr1', from => '1', to => '100000' }, { chrom => 'fake_chr1', from => '100001', to => '200000' }, { chrom => 'fake_chr1', from => '200001', to => '290640' }, { chrom => 'fake_chr2', from => '1', to => '100000' }, { chrom => 'fake_chr2', from => '100001', to => '200000' }, { chrom => 'fake_chr2', from => '200001', to => '300000' }, { chrom => 'fake_chr2', from => '300001', to => '400000' }, { chrom => 'fake_chr2', from => '400001', to => '500000' }, { chrom => 'fake_chr2', from => '500001', to => '600000' }, { chrom => 'fake_chr2', from => '600001', to => '700000' }, { chrom => 'fake_chr2', from => '700001', to => '800000' }, { chrom => 'fake_chr2', from => '800001', to => '900000' }, { chrom => 'fake_chr2', from => '900001', to => '1000000' }, { chrom => 'fake_chr2', from => '1000001', to => '1100000' }, { chrom => 'fake_chr2', from => '1100001', to => '1200000' }, { chrom => 'fake_chr2', from => '1200001', to => '1300000' }, { chrom => 'fake_chr2', from => '1300001', to => '1400000' }, { chrom => 'fake_chr2', from => '1400001', to => '1500000' }, { chrom => 'fake_chr2', from => '1500001', to => '1600000' }, { chrom => 'fake_chr2', from => '1600001', to => '1700000' }, { chrom => 'fake_chr2', from => '1700001', to => '1716851' }];

my @calling_files;
my $element_id = 1;
foreach my $sample (qw(1 2)) {
    foreach my $chunk (@$chunks) {
        my @output_subdirs = output_subdirs($element_id, 1);
        my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
        push(@calling_files, file(@output_subdirs, '4_gatk_unified_genotyper', qq[$region.gatk.vcf.gz]));
        push(@calling_files, file(@output_subdirs, '4_gatk_unified_genotyper', qq[$region.gatk.vcf.gz.tbi]));
        ++$element_id;
    }
}

ok handle_pipeline(), 'calling and concat pipelines ran ok';

my @concat_files;
foreach my $sample (qw(1 2)) {
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@concat_files, file(@output_subdirs, '1_vcf_concat', qq[merged.vcf.gz]));
    push(@concat_files, file(@output_subdirs, '1_vcf_concat', qq[merged.vcf.gz.tbi]));
    ++$element_id;
}

ok handle_pipeline(@calling_files, @concat_files), 'all expected calling and concat files were created';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 4, dataelement => 1)->cmd_summary->summary], ['java $jvm_args -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference_fasta -I $bams_list -o $vcf_file --genotype_likelihoods_model BOTH -stand_call_conf 0.0 -stand_emit_conf 0.0 -L $region'], 'cmd summaries for the major steps were as expected';

# check final vcfs have metadata
is_deeply [VRPipe::File->get(path => $concat_files[0])->metadata, VRPipe::File->get(path => $concat_files[2])->metadata],
  [{
        sample              => 'SAMPLE01',
        chunk_override_file => $override_file,
        caller              => 'GATK_UnifiedGenotyper'
    },
    {
        sample              => 'SAMPLE02',
        chunk_override_file => $override_file,
        caller              => 'GATK_UnifiedGenotyper'
    }
  ],
  'final merged vcfs have correct metadata';

finish;

