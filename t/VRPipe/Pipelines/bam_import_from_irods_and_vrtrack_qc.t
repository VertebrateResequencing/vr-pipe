#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 10;
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

# setup pipeline
my $output_dir = get_output_dir('bam_import_from_irods_and_vrtrack_qc');
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

ok my $import_qc_pipeline = VRPipe::Pipeline->get(name => 'bam_import_from_irods_and_vrtrack_qc'), 'able to get the bam_import_from_irods_and_vrtrack_qc pipeline';
my @s_names;
foreach my $stepmember ($import_qc_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename fasta_gc_stats bamcheck plot_bamcheck vrtrack_update_mapstats)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data pombe_ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$import_qc_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'pombe_ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

VRPipe::PipelineSetup->get(name => 'pombe import and qc',
                           datasource => $ds,
                           output_root => $output_dir,
                           pipeline => $import_qc_pipeline,
                           options => {reference_fasta => $ref_fa,
                                       reference_assembly_name => 'SPombe1',
                                       reference_public_url => 'ftp://s.pombe.com/ref.fa',
                                       reference_species => 'S.Pombe',
                                       bamcheck_options => '-q 20',
                                       exome_targets_file => file(qw(t data pombe_ref.fa.targets))->absolute->stringify,
                                       vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB},
                                       cleanup => 1});

my @irods_files;
my @qc_files;
my $element_id = 29;
foreach my $num (qw(1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31)) {
    my @output_subdirs = output_subdirs($element_id--);
    my $basename = '7369_5#'.$num;
    
    push(@irods_files, file($irods_dir, $basename.'.bam'));
    
    push(@qc_files, file(@output_subdirs, '3_bamcheck', $basename.'.bam.bamcheck'));
    foreach my $kind (qw(quals-hm quals quals2 quals3 insert-size gc-content gc-depth acgt-cycles coverage mism-per-cycle indel-dist indel-cycles)) {
        push(@qc_files, file(@output_subdirs, '4_plot_bamcheck', $basename.'-'.$kind.'.png'));
    }
}
ok handle_pipeline(@irods_files, @qc_files), 'irods import and qc graphs pipeline ran ok';

my $meta = VRPipe::File->get(path => $irods_files[0])->metadata;
is_deeply $meta, {bases => "0",
center_name => "SC",
expected_md5 => "8194f6c33299784d78e8d16fb05eb1c6",
individual => "SC_MFY5249218",
insert_size => "317",
lane => "7369_5#1",
lane_id => "85",
library => "4103711",
mate => "",
paired => "1",
platform => "SLX",
population => "Population",
reads => "0",
sample => "SC_MFY5249218",
sample_id => "",
study => "ERP001017",
study_name => "ERP001017",
targeted_avg_read_length => "100",
targeted_bases => "780614400",
targeted_bases_mapped => "778639100",
targeted_bases_mapped_c => "777515960",
targeted_bases_of_100X_coverage => "170588",
targeted_bases_of_10X_coverage => "11381918",
targeted_bases_of_1X_coverage => "11420817",
targeted_bases_of_20X_coverage => "11232925",
targeted_bases_of_2X_coverage => "11419827",
targeted_bases_of_50X_coverage => "9423004",
targeted_bases_of_5X_coverage => "11413494",
targeted_bases_trimmed => "1743015",
targeted_error_rate => "3.555864e-03",
targeted_forward_reads => "3903092",
targeted_mean_coverage => "65.47",
targeted_mean_insert_size => "313.3",
targeted_paired => "1",
targeted_reads => "7806144",
targeted_reads_mapped => "7786391",
targeted_reads_paired => "7763595",
targeted_reverse_reads => "3903052",
targeted_rmdup_bases => "772612200",
targeted_rmdup_bases_mapped => "770636900",
targeted_rmdup_reads => "7726122",
targeted_rmdup_reads_mapped => "7706369",
targeted_sd_insert_size => "80.8",
withdrawn => "0"}, 'metadata correct for one of the bam files';


# we'll also test the whole chain of pipelines we typically run in
# VertebrateResequencing at the Sanger: improvement followed by genotype check
# & auto_qc followed by merge up to the sample level

# improve
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
                                    datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                          method => 'all',
                                                                          source => 'pombe import and qc[1:local_files]',
                                                                          options => { }),
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

handle_pipeline();

my @improved_bams;
$element_id = 58;
foreach my $num (qw(1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31)) {
    my @output_subdirs = output_subdirs($element_id--, 2);
    push(@improved_bams, file(@output_subdirs, '10_bam_reheader', '7369_5#'.$num.'.realign.recal.calmd.bam'));
}
ok handle_pipeline(@improved_bams), 'chained improvement pipeline ran ok';

$meta = VRPipe::File->get(path => $improved_bams[0])->metadata;
my $opc = delete $meta->{original_pg_chain};
is_deeply $meta, {avg_read_length => "100",
bases => "1000298000",
bases_mapped => "990021900",
bases_mapped_c => "988611944",
bases_of_100X_coverage => "241872",
bases_of_10X_coverage => "12566325",
bases_of_1X_coverage => "12608642",
bases_of_20X_coverage => "12402778",
bases_of_2X_coverage => "12607594",
bases_of_50X_coverage => "10417019",
bases_of_5X_coverage => "12600810",
bases_trimmed => "0",
center_name => "SC",
error_rate => "3.613777e-03",
expected_md5 => "8194f6c33299784d78e8d16fb05eb1c6",
forward_reads => "5001490",
individual => "SC_MFY5249218",
insert_size => "317",
lane => "7369_5#1",
lane_id => "85",
library => "4103711",
mean_coverage => "65.99",
mean_insert_size => "311.8",
paired => "1",
platform => "SLX",
population => "Population",
reads => "10002980",
reads_mapped => "9900219",
reads_paired => "9867534",
reverse_reads => "5001490",
rmdup_bases => "982264700",
rmdup_bases_mapped => "971988600",
rmdup_reads => "9822647",
rmdup_reads_mapped => "9719886",
sample => "SC_MFY5249218",
sd_insert_size => "81.0",
study => "ERP001017",
study_name => "ERP001017",
targeted_avg_read_length => "100",
targeted_bases => "780614400",
targeted_bases_mapped => "778639100",
targeted_bases_mapped_c => "777515960",
targeted_bases_of_100X_coverage => "170588",
targeted_bases_of_10X_coverage => "11381918",
targeted_bases_of_1X_coverage => "11420817",
targeted_bases_of_20X_coverage => "11232925",
targeted_bases_of_2X_coverage => "11419827",
targeted_bases_of_50X_coverage => "9423004",
targeted_bases_of_5X_coverage => "11413494",
targeted_bases_trimmed => "1743015",
targeted_error_rate => "3.555864e-03",
targeted_forward_reads => "3903092",
targeted_mean_coverage => "65.47",
targeted_mean_insert_size => "313.3",
targeted_paired => "1",
targeted_reads => "7806144",
targeted_reads_mapped => "7786391",
targeted_reads_paired => "7763595",
targeted_reverse_reads => "3903052",
targeted_rmdup_bases => "772612200",
targeted_rmdup_bases_mapped => "770636900",
targeted_rmdup_reads => "7726122",
targeted_rmdup_reads_mapped => "7706369",
targeted_sd_insert_size => "80.8",
withdrawn => "0"}, 'metadata correct for one of the improved bam files: '.VRPipe::File->get(path => $improved_bams[0])->path;
$meta->{original_pg_chain} = $opc;


#*** genotype check that writes results to vrtrack pipeline...


# mergeup
my $merge_lanes_pipelinesetup = VRPipe::PipelineSetup->get(name => 'pombe merge lanes',
                                                           datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                                 method => 'group_by_metadata',
                                                                                                 source => 'pombe improvement[10]',
                                                                                                 options => { metadata_keys => 'population|sample|platform|library' } ),
                                                           output_root => $output_dir,
                                                           pipeline => VRPipe::Pipeline->get(name => 'merge_lanes'),
                                                           options => { bam_tags_to_strip => 'OQ XM XG XO',
                                                                        bam_merge_keep_single_paired_separate => 1,
                                                                        bam_merge_memory => 200,
                                                                        cleanup => 1 });

my $merge_libraries_pipelinesetup = VRPipe::PipelineSetup->get(name => 'pombe merge libraries',
                                                               datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                                     method => 'group_by_metadata',
                                                                                                     source => 'pombe merge lanes[3:markdup_bam_files]',
                                                                                                     options => { metadata_keys => 'population|sample|platform' } ),
                                                               output_root => $output_dir,
                                                               pipeline => VRPipe::Pipeline->get(name => 'merge_libraries_and_split'),
                                                               options => { bam_merge_keep_single_paired_separate => 0,
                                                                            reference_fasta => $ref_fa,
                                                                            reference_assembly_name => 'SPombe1',
                                                                            reference_public_url => 'ftp://s.pombe.com/ref.fa',
                                                                            reference_species => 'S.Pombe',
                                                                            bam_merge_memory => 200,
                                                                            split_bam_make_unmapped => 1,
                                                                            cleanup => 1,
                                                                            cleanup_inputs => 1,
                                                                            remove_merged_bams => 1 });

handle_pipeline();
my @merged_bams;
foreach my $element_id (89..117) {
    my @output_subdirs = output_subdirs($element_id, 4);
    foreach my $chrom (qw(chromIII chromI chromII chromAB325691 chromMT unmapped)) {
        push(@merged_bams, file(@output_subdirs, '4_bam_split_by_sequence', $chrom.'.pe.bam'));
    }
}
ok handle_pipeline(@merged_bams), 'chained mergeup pipeline ran ok';

finish;