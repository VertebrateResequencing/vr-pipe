#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 15;
    use VRPipeTest (
        required_env    => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe    => [qw(samtools bcftools)],
        required_module => [qw(Vcf)]
    );
    use TestPipelines;
}

ok my $calling_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_mpileup_via_bcf'), 'able to get the snp_calling_mpileup_via_bcf pipeline';

my $calling_dir = get_output_dir('wgs_multi_sample_mpileup_calling_test');

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
    pipeline    => $calling_pipeline,
    options     => {
        reference_fasta          => $ref_fa,
        samtools_mpileup_options => '-EDSV -C50 -m2 -F0.0005 -d 10000 -g',
        cleanup                  => 0,
    }
);

my $chunks = [{ chrom => 11, from => 1, to => 10000000 }, { chrom => 11, from => 10000001, to => 20000000 }, { chrom => 11, from => 20000001, to => 30000000 }, { chrom => 11, from => 30000001, to => 40000000 }, { chrom => 11, from => 40000001, to => 50000000 }, { chrom => 11, from => 50000001, to => 60000000 }, { chrom => 11, from => 60000001, to => 70000000 }, { chrom => 11, from => 70000001, to => 80000000 }, { chrom => 11, from => 80000001, to => 90000000 }, { chrom => 11, from => 90000001, to => 100000000 }, { chrom => 11, from => 100000001, to => 110000000 }, { chrom => 11, from => 110000001, to => 120000000 }, { chrom => 11, from => 120000001, to => 130000000 }, { chrom => 11, from => 130000001, to => 135006516 }, { chrom => 20, from => 1, to => 10000000 }, { chrom => 20, from => 10000001, to => 20000000 }, { chrom => 20, from => 20000001, to => 30000000 }, { chrom => 20, from => 30000001, to => 40000000 }, { chrom => 20, from => 40000001, to => 50000000 }, { chrom => 20, from => 50000001, to => 60000000 }, { chrom => 20, from => 60000001, to => 63025520 }];

my $element_id = 1;
my @calling_files;
my @restart_files;
foreach my $chunk (@$chunks) {
    my @output_subdirs = output_subdirs($element_id, 1);
    my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
    push(@calling_files, file(@output_subdirs, '2_mpileup_bcf', qq[$region.mpileup.bcf]));
    push(@calling_files, file(@output_subdirs, '3_bcf_to_vcf',  qq[$region.mpileup.vcf.gz]));
    push(@calling_files, file(@output_subdirs, '3_bcf_to_vcf',  qq[$region.mpileup.vcf.gz.tbi]));
    push(@restart_files, ($calling_files[-1], $calling_files[-2], $calling_files[-3])) if ($region eq '20_1-10000000');
    ++$element_id;
}

ok handle_pipeline(@calling_files), 'pipeline ran ok and all calling files were created';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 2, dataelement => 1)->cmd_summary->summary, VRPipe::StepState->get(pipelinesetup => 1, stepmember => 3, dataelement => 1)->cmd_summary->summary], ['samtools mpileup -EDSV -C50 -m2 -F0.0005 -d 10000 -g -r $region -f $reference_fasta -b $bams_list > $bcf_file', 'bcftools view -p 0.99 -vcgN -s $samples_file $bcf_file | bgzip -c > $vcf_file'], 'cmd summaries for the major steps were as expected';

# check override options work
# indel calling is turned off for chromosome 20 by options in override file
# check that no indels called on chrom20
cmp_ok indel_count($restart_files[1]->absolute->stringify), '>', 0, 'empty chunk_override_file worked as expected';
$original_override_file->copy($override_file);
VRPipe::DataElementState->get(pipelinesetup => 1, dataelement => 15)->start_from_scratch();
my @restart_file_exist = map { -s $_->absolute->stringify ? 1 : 0 } @restart_files;
is_deeply \@restart_file_exist, [0, 0, 0], 'correct files were removed after start from scratch';
ok handle_pipeline(@calling_files), 'pipeline ran ok after element started from scratch';
is indel_count($restart_files[1]->absolute->stringify), 0, 'chunk_override_file correctly overrode chunk option';

# run annotation pipeline
my $vep_cache = file(qw(t data vep_cache))->absolute->stringify;

ok my $annotation_pipeline = VRPipe::Pipeline->create(name => 'vcf_vep_annotate'), 'able to get the vcf_vep_annotate pipeline';
ok my $concat_pipeline     = VRPipe::Pipeline->create(name => 'vcf_concat'),       'able to get the vcf_concat pipeline';

VRPipe::PipelineSetup->create(
    name       => 'wgs test annotation',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => '1[3]',
        options => {}
    ),
    output_root => $calling_dir,
    pipeline    => $annotation_pipeline,
    options     => {
        vep_options                => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir $vep_cache",
        'vcf2consequences_options' => "--grantham",
        cleanup                    => 0
    }
);

VRPipe::PipelineSetup->create(
    name       => 'concat vcfs',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '2[2]',
        options => { metadata_keys => 'analysis_group|chrom' },
    ),
    output_root => $calling_dir,
    pipeline    => $concat_pipeline,
    options     => { cleanup => 0 }
);

ok handle_pipeline(), 'annotation and concat pipelines ran okay';

my @annotation_files;
foreach my $chunk (@$chunks) {
    my @output_subdirs = output_subdirs($element_id, 2);
    my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
    push(@annotation_files, file(@output_subdirs, '1_vep_analysis',         qq[$region.mpileup.vep.txt]));
    push(@annotation_files, file(@output_subdirs, '2_vcf_vep_consequences', qq[$region.mpileup.conseq.vcf.gz]));
    push(@annotation_files, file(@output_subdirs, '2_vcf_vep_consequences', qq[$region.mpileup.conseq.vcf.gz.tbi]));
    ++$element_id;
}

my @final_files;
foreach my $chrom (qw(11 20)) {
    my @output_subdirs = output_subdirs($element_id, 3);
    push(@final_files, file(@output_subdirs, '1_vcf_concat', 'merged.vcf.gz'));
    push(@final_files, file(@output_subdirs, '1_vcf_concat', 'merged.vcf.gz.tbi'));
    ++$element_id;
}

ok handle_pipeline(@annotation_files, @final_files), 'annotation and concat pipelines created all expected output files';

# check final vcfs have metadata
is_deeply [VRPipe::File->get(path => $final_files[0])->metadata, VRPipe::File->get(path => $final_files[2])->metadata],
  [{
        chrom               => '11',
        analysis_group      => 'low_coverage',
        chunk_override_file => $override_file->path->stringify,
        caller              => 'samtools_mpileup_bcftools'
    },
    {
        chrom               => '20',
        analysis_group      => 'low_coverage',
        chunk_override_file => $override_file->path->stringify,
        caller              => 'samtools_mpileup_bcftools'
    }
  ],
  'final merged vcfs have correct metadata';

# Call on subsets - EUR, ASN
my $ceu_samples_file = file(qw(t data wgs_calling_ceu.samples))->absolute->stringify;
my $asn_samples_file = file(qw(t data wgs_calling_asn.samples))->absolute->stringify;

ok my $bcf_calling_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_mpileup_from_bcf'), 'able to get the snp_calling_mpileup_from_bcf pipeline';

VRPipe::PipelineSetup->create(
    name       => 'wgs test CEU calling',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => '1[2]',  # bcf files
        options => {},
    ),
    output_root => $calling_dir,
    pipeline    => $bcf_calling_pipeline,
    options     => {
        reference_fasta => $ref_fa,
        sample_sex_file => $ceu_samples_file,
        cleanup         => 0,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'wgs test ASN calling',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => '1[2]',  # bcf files
        options => {},
    ),
    output_root => $calling_dir,
    pipeline    => $bcf_calling_pipeline,
    options     => {
        reference_fasta => $ref_fa,
        sample_sex_file => $asn_samples_file,
        cleanup         => 0,
    }
);

# Recall on merged site list
ok my $merge_pipeline = VRPipe::Pipeline->create(name => 'merge_vcfs_to_site_list_and_recall_from_bcf'), 'able to get the merge_vcfs_to_site_list pipeline';

VRPipe::PipelineSetup->create(
    name       => 'merge continental vcfs and recall from merged list',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '1[2]|4[1]|5[1]',                                  # bcf files, continental vcf files
        options => { metadata_keys => 'analysis_group|chrom|from|to' }
    ),
    output_root => $calling_dir,
    pipeline    => $merge_pipeline,
    options     => {
        reference_fasta => $ref_fa,
        sample_sex_file => $all_samples_file,
        cleanup         => 0
    }
);

VRPipe::PipelineSetup->create(
    name       => 'concat vcfs',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => '6[2]',
        options => { metadata_keys => 'analysis_group|chrom' },
    ),
    output_root => $calling_dir,
    pipeline    => $concat_pipeline,
    options     => { cleanup => 0 }
);

ok handle_pipeline(), 'continental calling and recalling on merged sites ran okay';

sub indel_count {
    my $vcf_path    = shift;
    my $indel_count = 0;
    my $vcf         = Vcf->new(file => $vcf_path);
    $vcf->parse_header();
    while (my $x = $vcf->next_data_array()) {
        my @alt = split(/,/, $$x[4]);
        for my $allele (@alt) {
            my ($type, $len) = $vcf->event_type({ REF => $$x[3] }, $allele);
            if ($type eq 'i') { ++$indel_count; }
        }
    }
    $vcf->close;
    return $indel_count;
}

finish;
