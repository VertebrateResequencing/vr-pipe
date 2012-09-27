#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK)],
        required_exe => [qw(vcf-concat samtools tabix)]
    );
    use TestPipelines;
}

ok my $calling_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_gatk_haplotype_caller'), 'able to get the snp_calling_gatk_haplotype_caller pipeline';
ok my $concat_pipeline = VRPipe::Pipeline->create(name => 'vcf_concat'), 'able to get the vcf_concat pipeline';

my $calling_dir = get_output_dir('haplotype_caller_test');

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

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data calling_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($calling_dir, 'calling_datasource.fofn'));
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
    name       => 'haplotype caller test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata_with_genome_chunking',
        method  => 'group_all',
        source  => $fofn_file->path,
        options => {
            reference_index     => file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify,
            chrom_list          => '11 20',
            chunk_size          => 10000000,
            chunk_overlap       => 0,
            chunk_override_file => $override_file,
        },
    ),
    output_root => $calling_dir,
    pipeline    => $calling_pipeline,
    options     => {
        reference_assembly_name  => 'NCBI37',
        reference_species        => 'Human',
        reference_fasta          => $ref_fa,
        haplotype_caller_options => '--genotyping_mode DISCOVERY -stand_call_conf 0.0 -stand_emit_conf 0.0',
        cleanup                  => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'vcf-concat test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '1[4]',
        options => { metadata_keys => 'chrom' }
    ),
    output_root => $calling_dir,
    pipeline    => $concat_pipeline,
    options     => {}
);

my $chunks = [{ chrom => 11, from => 1, to => 10000000 }, { chrom => 11, from => 10000001, to => 20000000 }, { chrom => 11, from => 20000001, to => 30000000 }, { chrom => 11, from => 30000001, to => 40000000 }, { chrom => 11, from => 40000001, to => 50000000 }, { chrom => 11, from => 50000001, to => 60000000 }, { chrom => 11, from => 60000001, to => 70000000 }, { chrom => 11, from => 70000001, to => 80000000 }, { chrom => 11, from => 80000001, to => 90000000 }, { chrom => 11, from => 90000001, to => 100000000 }, { chrom => 11, from => 100000001, to => 110000000 }, { chrom => 11, from => 110000001, to => 120000000 }, { chrom => 11, from => 120000001, to => 130000000 }, { chrom => 11, from => 130000001, to => 135006516 }, { chrom => 20, from => 1, to => 10000000 }, { chrom => 20, from => 10000001, to => 20000000 }, { chrom => 20, from => 20000001, to => 30000000 }, { chrom => 20, from => 30000001, to => 40000000 }, { chrom => 20, from => 40000001, to => 50000000 }, { chrom => 20, from => 50000001, to => 60000000 }, { chrom => 20, from => 60000001, to => 63025520 }];

ok handle_pipeline(), 'pipeline ran ok';

my $element_id = 1;
my @calling_files;
foreach my $chunk (@$chunks) {
    my @output_subdirs = output_subdirs($element_id, 1);
    my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
    push(@calling_files, file(@output_subdirs, '4_gatk_haplotype_caller', qq[$region.gatk_haplotype.vcf.gz]));
    push(@calling_files, file(@output_subdirs, '4_gatk_haplotype_caller', qq[$region.gatk_haplotype.vcf.gz.tbi]));
    ++$element_id;
}

my @concat_files;
foreach my $chrom (qw(11 20)) {
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@concat_files, file(@output_subdirs, '1_vcf_concat', qq[merged.vcf.gz]));
    push(@concat_files, file(@output_subdirs, '1_vcf_concat', qq[merged.vcf.gz.tbi]));
    ++$element_id;
}

ok handle_pipeline(@calling_files, @concat_files), 'all output files were created';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 4, dataelement => 1)->cmd_summary->summary], ['java $jvm_args -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference_fasta -I $bams_list -o $vcf_file --genotyping_mode DISCOVERY -stand_call_conf 0.0 -stand_emit_conf 0.0 -L $region'], 'cmd summaries for the major steps were as expected';

# check final vcfs have metadata
is_deeply [VRPipe::File->get(path => $concat_files[0])->metadata, VRPipe::File->get(path => $concat_files[2])->metadata],
  [{
        chrom               => '11',
        platform            => 'ILLUMINA',
        chunk_override_file => $override_file,
        caller              => 'GATK_HaplotypeCaller'
    },
    {
        chrom               => '20',
        platform            => 'ILLUMINA',
        chunk_override_file => $override_file,
        caller              => 'GATK_HaplotypeCaller'
    }
  ],
  'final merged vcfs have correct metadata';

done_testing;
exit;
