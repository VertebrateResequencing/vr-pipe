#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],);
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_split_and_vep_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->create(name => 'vcf_split_and_vep_annotate'), 'able to get the vcf_split_and_vep_annotate pipeline';
my @s_names;

foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(vcf_index vcf_split vep_analysis vcf_vep_consequences vcf_concat vcf_index);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $vep_cache  = file(qw(t data vep_cache))->absolute->stringify;
my $gerp_cache = file(qw(t data gerp_cache))->absolute->stringify;

# copy input vcfs to the output dir, since we will create .tbi files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data datasource.vcf_fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($output_dir, 'datasource.vcf_fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
while (<$ifh>) {
    chomp;
    my ($source_path, @meta) = split(/\t/, $_);
    my $source = file($source_path);
    my $dest = file($output_dir, $source->basename);
    copy($source, $dest);
    print $ofh join("\t", $dest, @meta);
    print $ofh "\n";
}
$orig_fofn_file->close;
$fofn_file->close;

my $test_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 'my vcf_chunked_vep_annotate pipeline setup',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => $fofn_file->path
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        cleanup                    => 0,
        reference_index            => file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify,
        chunk_size                 => 50000000,
        'vep_options'              => "--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir $vep_cache",
        'vcf2consequences_options' => "--grantham --gerp $gerp_cache",
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $file (qw(test1 test2)) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    foreach my $region (qw(1_1-50000000 1_50000001-100000000 1_100000001-150000000 1_150000001-200000000 1_200000001-249250621)) {
        push(@output_files, file(@output_subdirs, '2_vcf_split',            "$region.$file.vcf.gz"));
        push(@output_files, file(@output_subdirs, '3_vep_analysis',         "$region.$file.vep.txt"));
        push(@output_files, file(@output_subdirs, '4_vcf_vep_consequences', "$region.$file.conseq.vcf.gz"));
    }
    push(@final_files, file(@output_subdirs, '5_vcf_concat', 'merged.vcf.gz'));
    push(@final_files, file(@output_subdirs, '5_vcf_concat', 'merged.vcf.gz.tbi'));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

# check final vcfs have metadata
is_deeply [VRPipe::File->get(path => $final_files[0])->metadata, VRPipe::File->get(path => $final_files[2])->metadata], [{ chrom => '1' }, { chrom => '1' }], 'final merged vcfs have correct metadata';

done_testing;
exit;
