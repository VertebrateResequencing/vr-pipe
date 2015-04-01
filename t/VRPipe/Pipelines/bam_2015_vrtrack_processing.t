#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 5;
    use VRPipeTest;
    use VRPipeTest (
        # these are author-only tests (Sanger-specific)
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB VRPIPE_AUTHOR_TESTS VRPIPE_MOUSE_REF_DIR WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER SAMTOOLS GATK_JAVA)],
        required_exe => [qw(imeta)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

my $output_dir = get_output_dir('bam_2015_vrtrack_processing_no_recal');

# some mouse cram files
ok my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'uk10k',
    options => {
        file_query     => q[study_id = 2547 and type = cram and target = 1 and manual_qc like "%"],
        local_root_dir => $output_dir
    }
  ),
  'could create an irods datasource for cram files';

create_fresh_vrtrack_test_db();

ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_2015_vrtrack_processing_no_recal'), 'able to get the bam_2015_vrtrack_processing_no_recal pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names,
  [
    'vrtrack_populate_from_vrpipe_metadata',
    'irods_get_files_by_basename',
    'samtools_fasta_gc_stats',
    'samtools_bam_stats',
    'plot_bamstats',
    'vrtrack_update_mapstats',
    'sequence_dictionary',
    'bam_metadata',
    'bam_index',
    'gatk_target_interval_creator',
    'bam_realignment_around_known_indels',
    'bam_calculate_bq',
    'bam_reheader',
    'vrtrack_update_improved',
    'bam_index',
    'bcftools_generate_sites_file',
    'mpileup_vcf',
    'vcf_index',
    'bcftools_gtcheck',
    'bcftools_genotype_analysis',
    'vrtrack_auto_qc',
    'archive_files'
  ],
  'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data GRCm38.MT.ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $known_indels_source = file(qw(t data GRCm38.MT.known_indels.vcf));
my $known_indels = file($ref_dir, 'known_indels.vcf')->stringify;
copy($known_indels_source, $known_indels);

my $geno_source = file($ENV{VRPIPE_MOUSE_REF_DIR}, 'resources', 'genotypes', 'mouse_GRCm38_genotypes.vcf.gz');
my $geno = file($ref_dir, 'genotypes.vcf.gz')->stringify;
copy($geno_source,       $geno);
copy("$geno_source.tbi", "$geno.tbi");

my $dpf = VRPipe::File->create(path => file($output_dir, 'disc_pool_file'), type => 'txt');
my $dpfh = $dpf->openw;
my @pools;
foreach my $pool (qw(pool1 pool2 pool3)) {
    my $dir = dir($output_dir, 'archive_pools', $pool);
    $pipeline->make_path($dir);
    print $dpfh $dir, "\n";
    push(@pools, $dir);
}
$dpf->close;
my $pool_regex = join('|', @pools);

VRPipe::PipelineSetup->create(
    name        => 'mouse processing',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db                    => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir            => $output_dir,
        reference_fasta               => $ref_fa,
        reference_assembly_name       => 'GRCm38',
        reference_species             => 'Mus musculus',
        samtools_stats_options        => '-q 20',
        irods_convert_cram_to_bam     => file($ENV{SAMTOOLS}, 'samtools')->absolute->stringify,
        known_indels_for_realignment  => "-known $known_indels",
        gatk_count_covariates_options => '-l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
        genotypes_file                => $geno,
        # vcf_sample_from_metadata          => 'sample:public_name+sample',
        # expected_sample_from_metadata_key => 'public_name+sample',
        disc_pool_file => $dpf->path->stringify,
        gatk_path      => $ENV{GATK},
        picard_path    => $ENV{PICARD},
        java_exe       => $ENV{GATK_JAVA},
        cleanup        => 1
    }
);

my @expected_output_files;
#...
ok handle_pipeline(@expected_output_files), 'bam_2015_vrtrack_processing_no_recal pipeline ran ok for mouse cram files';

#*** needs some proper tests to ensure it all worked properly, but this is
# sufficient to know that we improved a bam and got stats in the VRTrack db

finish;
exit;

sub create_fresh_vrtrack_test_db {
    my %args = @_;
    
    my %cd = VRTrack::Factory->connection_details('rw');
    open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
    print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    foreach my $sql (VRTrack::VRTrack->schema()) {
        print $mysqlfh $sql;
    }
    
    close($mysqlfh);
}
