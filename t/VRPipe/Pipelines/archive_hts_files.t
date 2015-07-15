#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;
use Digest::MD5;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES');
    use TestPipelines;
}

my $archive_output_dir = get_output_dir('archive_hts_pipeline');
ok my $bam_pipeline = VRPipe::Pipeline->create(name => 'archive_bam_files'), 'able to get the archive_bam_files pipeline';
ok my $vcf_pipeline = VRPipe::Pipeline->create(name => 'archive_vcf_files'), 'able to get the archive_vcf_files pipeline';

my @s_names;
foreach my $stepmember ($bam_pipeline->step_members, $vcf_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(archive_bam_files archive_bam_index_files archive_vcf_files archive_vcf_index_files);
is_deeply \@s_names, \@expected_step_names, 'the pipelines have the correct steps';

# copy input files to the output dir, since we will be moving them
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data cram.fofn_with_metadata))->absolute);
my $fofn_file = VRPipe::File->create(path => file($archive_output_dir, 'cram.fofn_with_metadata'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
print $ofh scalar <$ifh>;
while (<$ifh>) {
    chomp;
    my ($source_path, @meta) = split(/\t/, $_);
    my $source = file($source_path);
    my $dest = file($archive_output_dir, $source->basename);
    copy($source, $dest);
    print $ofh join("\t", $dest, @meta);
    print $ofh "\n";
    copy("$source.crai", "$dest.crai");
    print $ofh join("\t", "$dest.crai", @meta);
    print $ofh "\n";
}
$orig_fofn_file->close;
$fofn_file->close;

# now run our archive_files pipeline, first creating pool dirs and the pool file
my $dpf = VRPipe::File->create(path => file($archive_output_dir, 'disc_pool_file'), type => 'txt');
my $dpfh = $dpf->openw;
my @pools;
foreach my $pool (qw(pool1 pool2 pool3)) {
    my $dir = dir($archive_output_dir, $pool);
    $bam_pipeline->make_path($dir);
    print $dpfh $dir, "\n";
    push(@pools, $dir);
}
$dpf->close;
my $pool_regex = join('|', @pools);

my $test_output_dir = get_output_dir('test_pipeline');
VRPipe::PipelineSetup->create(
    name       => 'cram archive test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => $fofn_file->path->stringify,
        options => { metadata_keys => 'lane' }
    ),
    output_root => $test_output_dir,
    pipeline    => $bam_pipeline,
    options     => { disc_pool_file => $dpf->path->stringify }
);

ok handle_pipeline(), 'archive_bam_files pipeline ran';

finish;
exit;

sub archive_file_location {
    my $file     = VRPipe::File->create(path => shift);
    my $dmd5     = Digest::MD5->new();
    my $ifile_id = $file->id;
    $dmd5->add($ifile_id);
    my $md5 = $dmd5->hexdigest;
    my @chars = split("", $md5);
    return (@chars[0 .. 3], $ifile_id, $file->basename);
}
