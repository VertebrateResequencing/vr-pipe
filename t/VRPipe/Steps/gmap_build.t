use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(gmap_build)]
    );
    use TestPipelines;
    use_ok('VRPipe::Steps::gmap_build');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('gmap_build', 'fasta_files');
is_deeply [$step->id, $step->description], [1, 'Indexes a reference genome fasta file, making it suitable for use in subsequent GSNAP mapping'], 'gmap_build step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'gmap_build',
    datasource => VRPipe::DataSource->create(
        type   => 'fofn',
        method => 'all',
        source => file(qw(t data gmap_build.fofn))->absolute
    ),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => { iit_file => file(qw(t data sacCer3_fake_snps.iit))->absolute, gmap_build_kmer_size => 15, gmap_build_fasta_files => file(qw(t data sacCer3 chr1.fa.gz))->absolute }
      #options     => {  gmap_build_kmer_size => 15, gmap_build_fasta_files => file(qw(t data sacCer3 chr1.fa.gz))->absolute }
);

my @output_subdirs = output_subdirs(1);
my $outputfile_1 = file(@output_subdirs, '1_gmap_build', 'mygenome', 'mygenome.chromosome');
my @outputfiles;
push(@outputfiles, $outputfile_1);

ok handle_pipeline(@outputfiles), 'gmap_build step ran ok, generating the expected output file';
warn $outputfile_1;
my $testfilecontents   = file(qw(t data gmap_build_sacCer3_chrI.chromosome))->slurp;
my $outputfilecontents = $outputfile_1->slurp;
is($testfilecontents, $outputfilecontents, "output file contain expected data");
finish;

