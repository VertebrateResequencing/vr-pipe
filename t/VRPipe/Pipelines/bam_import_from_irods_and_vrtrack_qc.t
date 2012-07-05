#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 21;
    # this test is Sanger-specific, only the author needs to run it
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
                    required_exe => [qw(iget iquest)]);
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

# setup a little VRTrack db that has its files in irods
my %cd = VRTrack::Factory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRPipe::File->get(path => file(qw(t data vrtrack_cerevisiae_wgs.sql))->absolute)->slurp;
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
is $results, 5, 'got correct number of bams from the vrtrack db';

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
my $element_id = 5;
my @lane_nums = qw(27 28 29 30 31);
foreach my $num (@lane_nums) {
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
withdrawn => "0",
population => "Population",
targeted_bases_of_20X_coverage => "11209858",
targeted_paired => "1",
targeted_reads_paired => "7819448",
targeted_bases_trimmed => "1766592",
individual => "SC_MFY5249244",
targeted_reverse_reads => "3930943",
targeted_bases_of_50X_coverage => "9399756",
sample => "SC_MFY5249244",
study => "ERP001017",
lane => "7369_5#27",
targeted_reads => "7861907",
targeted_rmdup_bases_mapped => "776426400",
targeted_sd_insert_size => "80.9",
targeted_error_rate => "3.597674e-03",
insert_size => "320",
targeted_rmdup_reads => "7783736",
targeted_bases_of_2X_coverage => "11417382",
paired => "1",
reads => "0",
project => "SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants",
targeted_mean_insert_size => "317.3",
library => "4103684",
targeted_bases_mapped => "784243500",
targeted_bases_of_10X_coverage => "11374735",
lane_id => "13",
targeted_bases_mapped_c => "783066516",
targeted_reads_mapped => "7842435",
targeted_bases_of_5X_coverage => "11410160",
targeted_bases => "786190700",
center_name => "SC",
platform => "SLX",
expected_md5 => "76ee5f08e761fa5baf637a6ad383ad19",
targeted_mean_coverage => "66.18",
targeted_avg_read_length => "100",
targeted_bases_of_1X_coverage => "11418447",
species => "Pombe",
targeted_rmdup_reads_mapped => "7764264",
targeted_bases_of_100X_coverage => "240676",
targeted_rmdup_bases => "778373600",
targeted_forward_reads => "3930964"}, 'metadata correct for one of the bam files';

my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
my $lane = VRTrack::Lane->new_by_name($vrtrack, '7369_5#27');
my $mapstats = $lane->latest_mapping;
is_deeply [$lane->is_processed('import'), $lane->is_processed('mapped'), $lane->is_processed('qc'), $mapstats->raw_reads], [1, 1, 1, 7861907], 'VRTrack database was updated correctly';


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

$output_dir = get_output_dir('bam_improvement');
VRPipe::PipelineSetup->get(name => 'pombe improvement',
                                    datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                          method => 'all',
                                                                          source => 'pombe import and qc[1]',
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
                                                delete_input_bams => 1,
                                                vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB}});

handle_pipeline();

my @improved_bams;
$element_id = 10;
foreach my $num (qw(27 28 29 30 31)) {
    my @output_subdirs = output_subdirs($element_id--, 2);
    push(@improved_bams, file(@output_subdirs, '10_bam_reheader', '7369_5#'.$num.'.realign.recal.calmd.bam'));
}
ok handle_pipeline(@improved_bams), 'chained improvement pipeline ran ok';

$meta = VRPipe::File->get(path => $improved_bams[0])->metadata;
my $opc = delete $meta->{original_pg_chain};
is_deeply $meta, {bases => "1155792800",
withdrawn => "0",
targeted_bases_of_20X_coverage => "11209858",
targeted_paired => "1",
individual => "SC_MFY5249244",
study => "ERP001017",
lane => "7369_5#27",
targeted_reads => "7861907",
targeted_rmdup_bases_mapped => "776426400",
targeted_sd_insert_size => "80.9",
bases_trimmed => "0",
targeted_error_rate => "3.597674e-03",
reads_mapped => "11460868",
insert_size => "320",
bases_of_50X_coverage => "10404533",
targeted_rmdup_reads => "7783736",
targeted_bases_of_2X_coverage => "11417382",
paired => "1",
reads => "11557928",
bases_mapped_c => "1144432834",
reverse_reads => "5778964",
targeted_bases_mapped => "784243500",
targeted_bases_of_10X_coverage => "11374735",
targeted_bases_mapped_c => "783066516",
rmdup_bases => "1108554300",
targeted_reads_mapped => "7842435",
targeted_bases_of_5X_coverage => "11410160",
center_name => "SC",
platform => "SLX",
expected_md5 => "76ee5f08e761fa5baf637a6ad383ad19",
species => "Pombe",
targeted_bases_of_100X_coverage => "240676",
targeted_rmdup_bases => "778373600",
targeted_forward_reads => "3930964",
population => "Population",
targeted_reads_paired => "7819448",
sd_insert_size => "81.3",
targeted_bases_trimmed => "1766592",
mean_insert_size => "314.4",
mean_coverage => "66.94",
bases_of_5X_coverage => "12597313",
bases_of_100X_coverage => "334539",
targeted_reverse_reads => "3930943",
reads_paired => "11424028",
targeted_bases_of_50X_coverage => "9399756",
sample => "SC_MFY5249244",
rmdup_reads => "11085543",
bases_of_2X_coverage => "12605219",
rmdup_bases_mapped => "1098848300",
bases_of_10X_coverage => "12559368",
bases_of_20X_coverage => "12379834",
bases_of_1X_coverage => "12606379",
project => "SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants",
targeted_mean_insert_size => "317.3",
library => "4103684",
avg_read_length => "100",
lane_id => "13",
targeted_bases => "786190700",
bases_mapped => "1146086800",
rmdup_reads_mapped => "10988483",
error_rate => "3.720504e-03",
targeted_mean_coverage => "66.18",
targeted_avg_read_length => "100",
targeted_bases_of_1X_coverage => "11418447",
forward_reads => "5778964",
targeted_rmdup_reads_mapped => "7764264",}, 'metadata correct for one of the improved bam files';
$meta->{original_pg_chain} = $opc;

$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
$lane = VRTrack::Lane->new_by_name($vrtrack, '7369_5#27');
$mapstats = $lane->latest_mapping;
is_deeply [$lane->is_processed('import'), $lane->is_processed('mapped'), $lane->is_processed('qc'), $lane->is_processed('improved'), $lane->raw_reads, $mapstats->raw_reads], [1, 1, 1, 1, 11557928, 7861907], 'VRTrack database was updated correctly after improvement';

# genotype check
# (the pombe_snps.bin was made by:
# samtools mpileup -q 20 -uf $ref *.bam | bcftools view -cgvI - | grep -v "^#" > vcf # (working on the import bams that have individual accessions in the headers)
# perl -Mstrict -we 'my %ofhs; my @s = qw(ERS074387 ERS074388 ERS074389 ERS074390 ERS074391); my %c = (ERS074387 => "SC_MFY5249244", ERS074388 => "SC_MFY5249245", ERS074389 => "SC_MFY5249246", ERS074390 => "SC_MFY5249247", ERS074391 => "SC_MFY5249248"); open(my $fh, "vcf"); while (<$fh>) { my ($chr, $pos, undef, $ref, $alt, undef, undef, undef, undef, @results) = split; next if $alt =~ /,/; foreach my $i (0..4) { my $s = $s[$i]; my $r = $results[$i]; ($r) = split(":", $r); my $al; if ($r eq "1/1") { $al = "$alt$alt"; } elsif ($r eq "0/1") { $al = "$ref$alt"; } elsif ($r eq "0/0") { $al = "$ref$ref"; } else { die "bad r $r at pos $pos\n"; } my $ind = $c{$s}; unless (defined $ofhs{$ind}) { open(my $ofh, ">$ind"); $ofhs{$ind} = $ofh; } my $ofh = $ofhs{$ind}; print $ofh "$chr $pos $al\n"; } }'
# ls SC_* > SC.fofn
# hapmap2bin -f SC.fofn > pombe_snps.bin
# )
my $snp_bin_source = file(qw(t data pombe_snps.bin));
my $snp_bin = file($res_dir, 'pombe_snps.bin')->stringify;
copy($snp_bin_source, $snp_bin);

$output_dir = get_output_dir('genotype_check');
VRPipe::PipelineSetup->get(name => 'genotype_checking',
			   datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                 method => 'all',
                                                                 source => 'pombe improvement[10]',
                                                                 options => { }),
			   output_root => $output_dir,
			   pipeline => VRPipe::Pipeline->get(name => 'bam_genotype_checking'),
			   options => {reference_fasta => $ref_fa,
			               hapmap2bin_sample_genotypes_file => $snp_bin,
				       expected_sample_metadata_key => 'sample'});

ok handle_pipeline(), 'bam_genotype_checking pipeline ran';
my $gtype_analysis_correct = 0;
foreach my $bam (@improved_bams) {
    my $meta = VRPipe::File->get(path => $bam)->metadata;
    my $ind = $meta->{individual} || next;
    my $ga = $meta->{gtype_analysis} || next;
    if ($ga =~ /^status=confirmed expected=$ind found=$ind/) {
        $gtype_analysis_correct++;
    }
}
is $gtype_analysis_correct, 5, 'all the bams had the expected gtype_analysis metadata added to them by the bam_genotype_checking pipeline';

# autoqc that writes results to vrtrack
$output_dir = get_output_dir('auto_qc');
my $auto_qc_ps = VRPipe::PipelineSetup->get(name => 'auto_qc',
                                            datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                  method => 'group_by_metadata',
                                                                                  source => 'pombe import and qc[3]|pombe improvement[10]',
                                                                                  options => { metadata_keys => 'lane' }),
                                            output_root => $output_dir,
                                            pipeline => VRPipe::Pipeline->get(name => 'vrtrack_auto_qc'),
                                            options => { vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB} });

ok handle_pipeline(), 'vrtrack_auto_qc pipeline ran';

my %actual_auto_qc_files;
my $lni = 0;
foreach my $de (@{$auto_qc_ps->datasource->elements}) {
    my @output_subdirs = output_subdirs($de->id, $auto_qc_ps->id);
    my $basename = '7369_5_'.$lane_nums[$lni++];
    my $aqcfile = VRPipe::File->get(path => file(@output_subdirs, '1_vrtrack_auto_qc', $basename.'.auto_qc.txt'));
    if ($aqcfile->s) {
        $actual_auto_qc_files{$basename} = [$aqcfile->slurp];
    }
}

my $passed_auto_qc_lanes = 0;
my $failed_auto_qc_libs = 0;
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
foreach my $num (@lane_nums) {
    my $lane = VRTrack::Lane->new_by_name($vrtrack, '7369_5#'.$num);
    $passed_auto_qc_lanes++ if $lane->auto_qc_status eq 'passed';
    my $lib = VRTrack::Library->new($vrtrack, $lane->library_id);
    $failed_auto_qc_libs++ if $lib->auto_qc_status eq 'failed';
}
is_deeply [$passed_auto_qc_lanes, $failed_auto_qc_libs], [5, 5], 'auto qc pipeline set all the lanes in VRTrack to passed, and their libraries to failed';

my $lane_gtype_correct = 0;
$passed_auto_qc_lanes = 0;
foreach my $bam (@improved_bams) {
    my $meta = VRPipe::File->get(path => $bam)->metadata;
    $passed_auto_qc_lanes++ if $meta->{auto_qc_status} eq 'passed';
    my $lane = $meta->{lane} || next;
    my $ind = $meta->{individual};
    my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $lane) || next;
    my $mapstats = $vrlane->latest_mapping || next;
    if ($mapstats->genotype_expected eq $ind && $mapstats->genotype_found eq $ind && $vrlane->genotype_status eq 'confirmed') {
        $lane_gtype_correct++;
    }
}
is_deeply [$lane_gtype_correct, $passed_auto_qc_lanes], [5, 5], 'auto qc also set bam metadata auto_qc_status to passed and updated VRTrack gt_status and mapstats gt_* for all lanes';

# test that we can re-run QC
$output_dir = get_output_dir('re_qc');
my $reqc_ps = VRPipe::PipelineSetup->get(name => 'pombe reqc',
                                         datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                               method => 'all',
                                                                               source => 'pombe improvement[10]',
                                                                               options => { }),
                                         output_root => $output_dir,
                                         pipeline => VRPipe::Pipeline->get(name => 'vrtrack_qc_graphs_and_auto_qc'),
                                         options => {reference_fasta => $ref_fa,
                                                     reference_assembly_name => 'SPombe1',
                                                     reference_public_url => 'ftp://s.pombe.com/ref.fa',
                                                     reference_species => 'S.Pombe',
                                                     bamcheck_options => '-q 30',
                                                     exome_targets_file => file(qw(t data pombe_ref.fa.targets))->absolute->stringify,
                                                     vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB},
                                                     auto_qc_error_rate => 0.01, # *** was trying to pick some option that makes some fail and other pass, but don't know what will do that...
                                                     cleanup => 1});

ok handle_pipeline(), 'vrtrack_qc_graphs_and_auto_qc pipeline ran';

$meta = VRPipe::File->get(path => $improved_bams[0])->metadata;
delete $meta->{original_pg_chain};
is_deeply $meta, {bases => "1155792800",
withdrawn => "0",
targeted_bases_of_20X_coverage => "11209858",
targeted_paired => "1",
individual => "SC_MFY5249244",
study => "ERP001017",
lane => "7369_5#27",
targeted_reads => "7861907",
targeted_rmdup_bases_mapped => "776426400",
targeted_sd_insert_size => "80.9",
bases_trimmed => "0",
targeted_error_rate => "3.597674e-03",
reads_mapped => "11460868",
insert_size => "320",
bases_of_50X_coverage => "10404533",
targeted_rmdup_reads => "7783736",
targeted_bases_of_2X_coverage => "11417382",
paired => "1",
reads => "11557928",
bases_mapped_c => "1144432834",
reverse_reads => "5778964",
targeted_bases_mapped => "784243500",
targeted_bases_of_10X_coverage => "11374735",
targeted_bases_mapped_c => "783066516",
rmdup_bases => "1108554300",
targeted_reads_mapped => "7842435",
targeted_bases_of_5X_coverage => "11410160",
center_name => "SC",
platform => "SLX",
expected_md5 => "76ee5f08e761fa5baf637a6ad383ad19",
species => "Pombe",
targeted_bases_of_100X_coverage => "240676",
targeted_rmdup_bases => "778373600",
targeted_forward_reads => "3930964",
population => "Population",
targeted_reads_paired => "7819448",
sd_insert_size => "81.3",
mean_insert_size => "314.4",
mean_coverage => "66.94",
bases_of_5X_coverage => "12597313",
bases_of_100X_coverage => "334539",
targeted_reverse_reads => "3930943",
reads_paired => "11424028",
targeted_bases_of_50X_coverage => "9399756",
sample => "SC_MFY5249244",
rmdup_reads => "11085543",
bases_of_2X_coverage => "12605219",
rmdup_bases_mapped => "1098848300",
bases_of_10X_coverage => "12559368",
bases_of_20X_coverage => "12379834",
bases_of_1X_coverage => "12606379",
project => "SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants",
library => "4103684",
avg_read_length => "100",
lane_id => "13",
targeted_bases => "786190700",
bases_mapped => "1146086800",
rmdup_reads_mapped => "10988483",
error_rate => "3.720504e-03",
targeted_mean_coverage => "66.18",
targeted_avg_read_length => "100",
targeted_bases_of_1X_coverage => "11418447",
forward_reads => "5778964",
targeted_rmdup_reads_mapped => "7764264",
auto_qc_status => 'passed',
gtype_analysis => 'status=confirmed expected=SC_MFY5249244 found=SC_MFY5249244 ratio=1.538',
targeted_bases_trimmed => "30289444",
targeted_mean_insert_size => "316.3", }, 'metadata got updated for one of the improved bam files after redoing QC';

$passed_auto_qc_lanes = 0;
foreach my $bam (@improved_bams) {
    my $meta = VRPipe::File->get(path => $bam)->metadata;
    $passed_auto_qc_lanes++ if $meta->{auto_qc_status} eq 'passed';
}
is $passed_auto_qc_lanes, 5, 'after redoing QC, still have 5 lanes passed';

# mergeup
$output_dir = get_output_dir('lane_merge');
VRPipe::PipelineSetup->get(name => 'pombe merge lanes',
                           datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                 method => 'group_by_metadata',
                                                                 source => 'pombe improvement[10]',
                                                                 options => { metadata_keys => 'population|sample|platform|library' } ),
                           output_root => $output_dir,
                           pipeline => VRPipe::Pipeline->get(name => 'bam_merge_lanes_and_fix_rgs'),
                           options => { bam_tags_to_strip => 'OQ XM XG XO',
                                        bam_merge_keep_single_paired_separate => 1,
                                        cleanup => 1,
                                        delete_input_bams => 1 });

$output_dir = get_output_dir('library_merge');
VRPipe::PipelineSetup->get(name => 'pombe merge libraries',
                           datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                 method => 'group_by_metadata',
                                                                 source => 'pombe merge lanes[4]',
                                                                 options => { metadata_keys => 'population|sample|platform' } ),
                           output_root => $output_dir,
                           pipeline => VRPipe::Pipeline->get(name => 'bam_merge_and_split'),
                           options => { bam_merge_keep_single_paired_separate => 0,
                                        split_bam_make_unmapped => 0,
                                        delete_input_bams => 1,
                                        remove_merged_bams => 1 });

# mergeacross
$output_dir = get_output_dir('merge_across');
my $mergeacross_ps = VRPipe::PipelineSetup->get(name => 'pombe mergeacross',
                                                datasource => VRPipe::DataSource->get(type => 'vrpipe',
                                                                                      method => 'group_by_metadata',
                                                                                      source => 'pombe merge libraries[2]',
                                                                                      options => { metadata_keys => 'split_sequence' } ),
                                                output_root => $output_dir,
                                                pipeline => VRPipe::Pipeline->get(name => 'bam_merge'),
                                                options => { bam_merge_keep_single_paired_separate => 0,
                                                             bam_merge_maximum_files => 2,
                                                             delete_input_bams => 1 });

handle_pipeline();
my @merged_bams;
foreach my $element (@{$mergeacross_ps->datasource->elements}) {
    my @output_subdirs = output_subdirs($element->id, $mergeacross_ps->id);
    for my $chunk (1..3) {
        push(@merged_bams, file(@output_subdirs, '1_bam_merge', "pe.$chunk.bam"));
    }
}
ok handle_pipeline(@merged_bams), 'chained mergeup -> mergeacross pipelines ran ok';

my %seen_splits;
foreach my $mbam (@merged_bams) {
    my $meta = VRPipe::File->get(path => $mbam)->metadata;
    my $split = $meta->{split_sequence} || next;
    $seen_splits{$split}++;
}
is_deeply \%seen_splits, {chromMT => 3, chromIII => 3, chromII => 3, chromI => 3, chromAB325691 => 3 }, 'there were 3 mergeacross bams for each chromosome, and the metadata was correct';

finish;
