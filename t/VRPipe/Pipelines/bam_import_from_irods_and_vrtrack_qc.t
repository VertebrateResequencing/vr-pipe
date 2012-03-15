#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 7;
    # this test is Sanger-specific, only the author needs to run it
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
                    required_exe => [qw(iget iquest)]);
    use TestPipelines;
    
    use_ok('VertRes::Utils::VRTrackFactory');
}

# setup a little VRTrack db that has its files in irods
my %cd = VertRes::Utils::VRTrackFactory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRPipe::File->get(path => file(qw(t data vrtrack_yeast_irods.sql))->absolute)->slurp;
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

my $output_dir = get_output_dir('bam_import_from_irods_and_vrtrack_qc_wgs');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

ok my $ds = VRPipe::DataSource->get(type => 'vrtrack',
                                    method => 'lane_bams',
                                    source => $ENV{VRPIPE_VRTRACK_TESTDB},
                                    options => {local_root_dir => $irods_dir}), 'could create a vrtrack datasource';
my $results = 0;
foreach my $element (@{$ds->elements}) {
    $results++;
}
is $results, 29, 'got correct number of bams from the vrtrack db';

ok my $import_qc_pipeline = VRPipe::Pipeline->get(name => 'bam_import_from_irods_and_vrtrack_qc_wgs'), 'able to get the bam_import_from_irods_and_vrtrack_qc_wgs pipeline';
my @s_names;
foreach my $stepmember ($import_qc_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename bamcheck plot_bamcheck bamcheck_rmdup bamcheck_stats_output vrtrack_update_mapstats)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data pombe_ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$import_qc_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'pombe_ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $stats_fa_source = $ref_fa_source.'.stats';
my $ref_fa_stats = $ref_fa.'.stats';
copy($stats_fa_source, $ref_fa_stats);

VRPipe::PipelineSetup->get(name => 'pombe import and qc',
                           datasource => $ds,
                           output_root => $output_dir,
                           pipeline => $import_qc_pipeline,
                           options => {reference_fasta => $ref_fa,
                                       reference_assembly_name => 'SPombe1',
                                       reference_public_url => 'ftp://s.pombe.com/ref.fa',
                                       reference_species => 'S.Pombe',
                                       reference_fasta_stats => $ref_fa_stats,
                                       bamcheck_options => '-q 20 -r', # -q20 -d for bamcheck_rmdup_options?
                                       vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB},
                                       cleanup => 1});

my @irods_files;
my @qc_files;
my $element_id = 29;
foreach my $num (qw(1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31)) {
    my @output_subdirs = output_subdirs($element_id--);
    my $basename = '7369_5#'.$num;
    
    push(@irods_files, file($irods_dir, $basename.'.bam'));
    
    push(@qc_files, file(@output_subdirs, '2_bamcheck', $basename.'.bam.bamcheck'));
    foreach my $kind (qw(quals-hm quals quals2 quals3 insert-size gc-content gc-depth acgt-cycles coverage mism-per-cycle indel-dist indel-cycles)) {
        push(@qc_files, file(@output_subdirs, '3_plot_bamcheck', $basename.'-'.$kind.'.png'));
    }
    push(@qc_files, file(@output_subdirs, '5_bamcheck_stats_output', $basename.'.bam.detailed_stats'));
}
ok handle_pipeline(@irods_files, @qc_files), 'irods import and qc graphs pipeline ran ok';

# we'll also test the whole chain of pipelines we typically run in
# VertebrateResequencing at the Sanger: improvement followed by genotype check
# & auto_qc followed by merge up to the sample level

my $res_dir = dir($output_dir, 'resources');
$import_qc_pipeline->make_path($res_dir);
my $known_indels_source = file(qw(t data known_indels_pombe.vcf.gz));
my $known_indels = file($res_dir, 'known_indels.vcf.gz')->stringify;
copy($known_indels_source, $known_indels);
copy($known_indels_source.'.tbi', $known_indels.'.tbi');

my $known_sites_source = file(qw(t data known_sites_pombe.vcf.gz));
my $known_sites = file($res_dir, 'known_sites.vcf.gz')->stringify;
copy($known_sites_source, $known_sites);
copy($known_sites_source.'.tbi', $known_sites.'.tbi');

VRPipe::PipelineSetup->get(name => 'pombe improvement',
                                    datasource => $ds,
                                    output_root => $output_dir,
                                    pipeline => VRPipe::Pipeline->get(name => 'bam_improvement_and_update_vrtrack'),
                                    options => {reference_fasta => $ref_fa,
                                                reference_assembly_name => 'SPombe1',
                                                reference_public_url => 'ftp://s.pombe.com/ref.fa',
                                                reference_species => 'S.Pombe',
                                                known_indels_for_realignment => "-known $known_indels",
                                                known_sites_for_recalibration => "-knownSites $known_sites",
                                                gatk_count_covariates_options => '-l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
                                                gatk_path => $ENV{GATK},
                                                picard_path => $ENV{PICARD},
                                                cleanup => 1,
                                                vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB}});

ok handle_pipeline, 'chained improvement pipeline ran ok';

finish;