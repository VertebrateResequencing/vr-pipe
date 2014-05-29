#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 2;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(samtools bcftools tabix)]
    );
    use TestPipelines;
}

# this is a cut-down version of wgs_single_sample_calling.t to just focus on
# some quick samtools | bcftools testing

ok my $mpileup_pipeline = VRPipe::Pipeline->create(name => 'snp_calling_mpileup'), 'able to get the snp_calling_mpileup pipeline';
my $calling_dir = get_output_dir('snp_calling_mpileup');

my $original_ref_fa = VRPipe::File->create(path => file(qw(t data S_suis_P17.fa))->absolute);
my $ref_fa = VRPipe::File->create(path => file($calling_dir, 'S_suis_P17.fa')->absolute);
$original_ref_fa->copy($ref_fa);

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
    name       => 'mpileup calling',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'all',
        source  => $fofn_file->path,
        options => {}
    ),
    output_root => $calling_dir,
    pipeline    => $mpileup_pipeline,
    options     => {
        reference_fasta => $ref_fa->path,
        cleanup         => 0,
    }
);

my (@output_files, @final_files);
my $element_id = 0;
foreach my $sample ('SAMPLE01', 'SAMPLE02') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '2_mpileup_vcf', "mpileup.vcf.gz"));
    push(@output_files, file(@output_subdirs, '2_mpileup_vcf', "mpileup.vcf.gz.tbi"));
}

ok handle_pipeline(@output_files), 'calling pipeline ran ok and created all expected output files';

finish;

