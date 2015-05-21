#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 21;
    # this test is Sanger-specific, only the author needs to run it
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
        required_exe => [qw(iget iquest)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

# setup a little VRTrack db that has its files in irods. The sql used here will
# need updating when VRTrack schema is incremented; to do that, put an exit
# after the following block and run this test to create the db, then update
# the db: update-vrtrack-schema -o -s ../vr-codebase/sql/VRTrack_schema_[...],
# then mysqldump it
my %cd = VRTrack::Factory->connection_details('rw');
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
my @sql = VRPipe::File->create(path => file(qw(t data vrtrack_cerevisiae_wgs.sql))->absolute)->slurp;
foreach my $sql (@sql) {
    print $mysqlfh $sql;
}
close($mysqlfh);

# setup pipeline
my $output_dir = get_output_dir('bam_import_from_irods_and_vrtrack_qc');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

ok my $ds = VRPipe::DataSource->create(
    type    => 'vrtrack',
    method  => 'lane_bams',
    source  => $ENV{VRPIPE_VRTRACK_TESTDB},
    options => { local_root_dir => $irods_dir }
  ),
  'could create a vrtrack datasource';
my $results = 0;
foreach my $element (@{ get_elements($ds) }) {
    $results++;
}
is $results, 5, 'got correct number of bams from the vrtrack db';

ok my $import_qc_pipeline = VRPipe::Pipeline->create(name => 'bam_import_from_irods_and_vrtrack_qc'), 'able to get the bam_import_from_irods_and_vrtrack_qc pipeline';
my @s_names;
foreach my $stepmember ($import_qc_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(irods_get_files_by_basename fasta_gc_stats bamcheck plot_bamcheck vrtrack_update_mapstats)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data pombe_ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$import_qc_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'pombe_ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

VRPipe::PipelineSetup->create(
    name        => 'pombe import and qc',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $import_qc_pipeline,
    options     => {
        reference_fasta         => $ref_fa,
        reference_assembly_name => 'SPombe1',
        reference_public_url    => 'ftp://s.pombe.com/ref.fa',
        reference_species       => 'S.Pombe',
        bamcheck_options        => '-q 20',
        exome_targets_file      => file(qw(t data pombe_ref.fa.targets))->absolute->stringify,
        vrtrack_db              => $ENV{VRPIPE_VRTRACK_TESTDB},
        cleanup                 => 1
    }
);

my @irods_files;
my @qc_files;
my $element_id = 5;
my @lane_nums  = qw(27 28 29 30 31);
foreach my $num (@lane_nums) {
    my @output_subdirs = output_subdirs($element_id--);
    my $basename       = '7369_5#' . $num;
    
    push(@irods_files, file($irods_dir, $basename . '.bam'));
    
    push(@qc_files, file(@output_subdirs, '3_bamcheck', $basename . '.bam.bamcheck'));
    foreach my $kind (qw(quals-hm quals quals2 quals3 insert-size gc-content gc-depth acgt-cycles coverage mism-per-cycle indel-dist indel-cycles)) {
        push(@qc_files, file(@output_subdirs, '4_plot_bamcheck', $basename . '-' . $kind . '.png'));
    }
}
ok handle_pipeline(@irods_files, @qc_files), 'irods import and qc graphs pipeline ran ok';

my $meta = VRPipe::File->get(path => $irods_files[0])->metadata;
is_deeply $meta,
  {
    bases                           => 1155792800,
    center_name                     => 'SC',
    expected_md5                    => '18dc7d5820b901dfba36be1e41603a1b',
    individual                      => 'SC_MFY5249244',
    insert_size                     => 320,
    lane                            => '7369_5#27',
    lane_id                         => 13,
    library                         => 4103684,
    paired                          => 1,
    platform                        => 'SLX',
    population                      => 'Population',
    project                         => 'SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants',
    reads                           => 11557928,
    sample                          => 'SC_MFY5249244',
    species                         => 'Pombe',
    study                           => 'ERP001017',
    targeted_avg_read_length        => 100,
    targeted_bases                  => 785559800,
    targeted_bases_mapped           => 783669700,
    targeted_bases_mapped_c         => 782746015,
    targeted_bases_of_100X_coverage => 236968,
    targeted_bases_of_10X_coverage  => 11374070,
    targeted_bases_of_1X_coverage   => 11418355,
    targeted_bases_of_20X_coverage  => 11207157,
    targeted_bases_of_2X_coverage   => 11417323,
    targeted_bases_of_50X_coverage  => 9395673,
    targeted_bases_of_5X_coverage   => 11409861,
    targeted_bases_trimmed          => 1738372,
    targeted_error_rate             => '3.641800e-03',
    targeted_filtered_reads         => 7855598,
    targeted_forward_reads          => 3927828,
    targeted_mean_coverage          => 66.13,
    targeted_mean_insert_size       => 317.3,
    targeted_paired                 => 1,
    targeted_reads                  => 7855598,
    targeted_reads_mapped           => 7836697,
    targeted_reads_paired           => 7817796,
    targeted_reverse_reads          => 3927770,
    targeted_rmdup_bases            => 777710400,
    targeted_rmdup_bases_mapped     => 774905867,
    targeted_rmdup_reads            => 7777104,
    targeted_rmdup_reads_mapped     => 7758203,
    targeted_sd_insert_size         => 80.9,
    withdrawn                       => 0
  },
  'metadata correct for one of the bam files';

my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
my $lane = VRTrack::Lane->new_by_name($vrtrack, '7369_5#27');
my $mapstats = $lane->latest_mapping;
is_deeply [$lane->is_processed('import'), $lane->is_processed('mapped'), $lane->is_processed('qc'), $mapstats->raw_reads], [1, 1, 1, 7855598], 'VRTrack database was updated correctly';

# we'll also test the whole chain of pipelines we typically run in
# VertebrateResequencing at the Sanger: improvement followed by genotype check
# & auto_qc followed by merge up to the sample level

# improve
my $res_dir = dir($output_dir, 'resources');
$import_qc_pipeline->make_path($res_dir);
my $known_indels_source = file(qw(t data known_indels_pombe.vcf.gz));
my $known_indels = file($res_dir, 'known_indels.vcf.gz')->stringify;
copy($known_indels_source,          $known_indels);
copy($known_indels_source . '.tbi', $known_indels . '.tbi');

my $known_sites_source = file(qw(t data known_sites_pombe.vcf.gz));
my $known_sites = file($res_dir, 'known_sites.vcf.gz')->stringify;
copy($known_sites_source,          $known_sites);
copy($known_sites_source . '.tbi', $known_sites . '.tbi');

$output_dir = get_output_dir('bam_improvement');
VRPipe::PipelineSetup->create(
    name       => 'pombe improvement',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'pombe import and qc[1]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'bam_improvement_and_update_vrtrack'),
    options     => {
        reference_fasta               => $ref_fa,
        reference_assembly_name       => 'SPombe1',
        reference_public_url          => 'ftp://s.pombe.com/ref.fa',
        reference_species             => 'S.Pombe',
        known_indels_for_realignment  => "-known $known_indels",
        known_sites_for_recalibration => "-knownSites $known_sites",
        gatk_count_covariates_options => '-l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
        gatk_path                     => $ENV{GATK},
        picard_path                   => $ENV{PICARD},
        cleanup                       => 1,
        delete_input_bams             => 1,
        vrtrack_db                    => $ENV{VRPIPE_VRTRACK_TESTDB}
    }
);

handle_pipeline();

my @improved_bams;
$element_id = 10;
foreach my $num (qw(27 28 29 30 31)) {
    my @output_subdirs = output_subdirs($element_id--, 2);
    push(@improved_bams, file(@output_subdirs, '10_bam_reheader', '7369_5#' . $num . '.realign.recal.calmd.bam'));
}
ok handle_pipeline(@improved_bams), 'chained improvement pipeline ran ok';

$meta = VRPipe::File->get(path => $improved_bams[0])->metadata;
warn "checking metadata of $improved_bams[0]\n";
my $opc = delete $meta->{original_pg_chain};
is_deeply $meta,
  {
    avg_read_length                 => 100,
    bases                           => 1155792800,
    bases_mapped                    => 1146563200,
    bases_mapped_c                  => 1145243236,
    bases_of_100X_coverage          => 332396,
    bases_of_10X_coverage           => 12577952,
    bases_of_1X_coverage            => 12626444,
    bases_of_20X_coverage           => 12394676,
    bases_of_2X_coverage            => 12625336,
    bases_of_50X_coverage           => 10410466,
    bases_of_5X_coverage            => 12616903,
    bases_trimmed                   => 0,
    center_name                     => 'SC',
    error_rate                      => '3.762939e-03',
    expected_md5                    => '18dc7d5820b901dfba36be1e41603a1b',
    filtered_reads                  => 11557928,
    forward_reads                   => 5778964,
    individual                      => 'SC_MFY5249244',
    insert_size                     => 320,
    lane                            => '7369_5#27',
    lane_id                         => 13,
    library                         => 4103684,
    mean_coverage                   => 66.88,
    mean_insert_size                => 314.4,
    paired                          => 1,
    platform                        => 'SLX',
    population                      => 'Population',
    project                         => 'SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants',
    reads                           => 11557928,
    reads_mapped                    => 11465632,
    reads_paired                    => 11434358,
    reverse_reads                   => 5778964,
    rmdup_bases                     => 1108486300,
    rmdup_bases_mapped              => 1097991197,
    rmdup_reads                     => 11084863,
    rmdup_reads_mapped              => 10992567,
    sample                          => 'SC_MFY5249244',
    sd_insert_size                  => 81.2,
    species                         => 'Pombe',
    study                           => 'ERP001017',
    targeted_avg_read_length        => 100,
    targeted_bases                  => 785559800,
    targeted_bases_mapped           => 783669700,
    targeted_bases_mapped_c         => 782746015,
    targeted_bases_of_100X_coverage => 236968,
    targeted_bases_of_10X_coverage  => 11374070,
    targeted_bases_of_1X_coverage   => 11418355,
    targeted_bases_of_20X_coverage  => 11207157,
    targeted_bases_of_2X_coverage   => 11417323,
    targeted_bases_of_50X_coverage  => 9395673,
    targeted_bases_of_5X_coverage   => 11409861,
    targeted_bases_trimmed          => 1738372,
    targeted_error_rate             => '3.641800e-03',
    targeted_filtered_reads         => 7855598,
    targeted_forward_reads          => 3927828,
    targeted_mean_coverage          => 66.13,
    targeted_mean_insert_size       => 316.3,
    targeted_paired                 => 1,
    targeted_reads                  => 7855598,
    targeted_reads_mapped           => 7836697,
    targeted_reads_paired           => 7817796,
    targeted_reverse_reads          => 3927770,
    targeted_rmdup_bases            => 777710400,
    targeted_rmdup_bases_mapped     => 774905867,
    targeted_rmdup_reads            => 7777104,
    targeted_rmdup_reads_mapped     => 7758203,
    targeted_sd_insert_size         => 80.9,
    withdrawn                       => 0
  },
  'metadata correct for one of the improved bam files';
$meta->{original_pg_chain} = $opc;

$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
$lane = VRTrack::Lane->new_by_name($vrtrack, '7369_5#27');
$mapstats = $lane->latest_mapping;
is_deeply [$lane->is_processed('import'), $lane->is_processed('mapped'), $lane->is_processed('qc'), $lane->is_processed('improved'), $lane->raw_reads, $mapstats->raw_reads], [1, 1, 1, 1, 11557928, 7855598], 'VRTrack database was updated correctly after improvement';

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
VRPipe::PipelineSetup->create(
    name       => 'genotype_checking',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'pombe improvement[10]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'bam_genotype_checking'),
    options     => {
        reference_fasta                  => $ref_fa,
        hapmap2bin_sample_genotypes_file => $snp_bin,
        expected_sample_metadata_key     => 'sample'
    }
);

ok handle_pipeline(), 'bam_genotype_checking pipeline ran';
my $gtype_analysis_correct = 0;
foreach my $bam (@improved_bams) {
    my $meta = VRPipe::File->get(path => $bam)->metadata;
    my $ind = $meta->{individual}     || next;
    my $ga  = $meta->{gtype_analysis} || next;
    if ($ga =~ /^status=confirmed expected=$ind found=$ind/) {
        $gtype_analysis_correct++;
    }
}
#*** samtools v1+ result in 4 of these failing, need to investigate why and if
# gtypex_genotype_analysis step needs adjusting
is $gtype_analysis_correct, 5, 'all the bams had the expected gtype_analysis metadata added to them by the bam_genotype_checking pipeline';

# autoqc that writes results to vrtrack
$output_dir = get_output_dir('auto_qc');
my $auto_qc_ps = VRPipe::PipelineSetup->create(
    name       => 'auto_qc',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'pombe import and qc[3]|pombe improvement[10]',
        options => { metadata_keys => 'lane' }
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'vrtrack_auto_qc'),
    options     => { vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB} }
);

ok handle_pipeline(), 'vrtrack_auto_qc pipeline ran';

my %actual_auto_qc_files;
my $lni = 0;
foreach my $de (@{ get_elements($auto_qc_ps->datasource) }) {
    my @output_subdirs = output_subdirs($de->id, $auto_qc_ps->id);
    my $basename       = '7369_5_' . $lane_nums[$lni++];
    my $aqcfile        = VRPipe::File->create(path => file(@output_subdirs, '1_vrtrack_auto_qc', $basename . '.auto_qc.txt'));
    if ($aqcfile->s) {
        $actual_auto_qc_files{$basename} = [$aqcfile->slurp];
    }
}

my $passed_auto_qc_lanes = 0;
my $failed_auto_qc_libs  = 0;
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
foreach my $num (@lane_nums) {
    my $lane = VRTrack::Lane->new_by_name($vrtrack, '7369_5#' . $num);
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
    my $lane     = $meta->{lane} || next;
    my $ind      = $meta->{individual};
    my $vrlane   = VRTrack::Lane->new_by_name($vrtrack, $lane) || next;
    my $mapstats = $vrlane->latest_mapping || next;
    if ($mapstats->genotype_expected eq $ind && $mapstats->genotype_found eq $ind && $vrlane->genotype_status eq 'confirmed') {
        $lane_gtype_correct++;
    }
}
is_deeply [$lane_gtype_correct, $passed_auto_qc_lanes], [5, 5], 'auto qc also set bam metadata auto_qc_status to passed and updated VRTrack gt_status and mapstats gt_* for all lanes';

# test that we can re-run QC
$output_dir = get_output_dir('re_qc');
my $reqc_ps = VRPipe::PipelineSetup->create(
    name       => 'pombe reqc',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'pombe improvement[10]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'vrtrack_qc_graphs_and_auto_qc'),
    options     => {
        reference_fasta         => $ref_fa,
        reference_assembly_name => 'SPombe1',
        reference_public_url    => 'ftp://s.pombe.com/ref.fa',
        reference_species       => 'S.Pombe',
        bamcheck_options        => '-q 30',
        exome_targets_file      => file(qw(t data pombe_ref.fa.targets))->absolute->stringify,
        vrtrack_db              => $ENV{VRPIPE_VRTRACK_TESTDB},
        auto_qc_error_rate      => 0.01,                                                      # *** was trying to pick some option that makes some fail and other pass, but don't know what will do that...
        cleanup                 => 1
    }
);

ok handle_pipeline(), 'vrtrack_qc_graphs_and_auto_qc pipeline ran';

$meta = VRPipe::File->get(path => $improved_bams[0])->metadata;
delete $meta->{original_pg_chain};
is_deeply $meta,
  {
    auto_qc_status                  => 'passed',
    avg_read_length                 => 100,
    bases                           => 1155792800,
    bases_mapped                    => 1146563200,
    bases_mapped_c                  => 1145243236,
    bases_of_100X_coverage          => 332396,
    bases_of_10X_coverage           => 12577952,
    bases_of_1X_coverage            => 12626444,
    bases_of_20X_coverage           => 12394676,
    bases_of_2X_coverage            => 12625336,
    bases_of_50X_coverage           => 10410466,
    bases_of_5X_coverage            => 12616903,
    bases_trimmed                   => 0,
    center_name                     => 'SC',
    error_rate                      => '3.762939e-03',
    expected_md5                    => '18dc7d5820b901dfba36be1e41603a1b',
    filtered_reads                  => 11557928,
    forward_reads                   => 5778964,
    gtype_analysis                  => 'status=confirmed expected=SC_MFY5249244 found=SC_MFY5249244 ratio=1.530',
    individual                      => 'SC_MFY5249244',
    insert_size                     => 320,
    lane                            => '7369_5#27',
    lane_id                         => 13,
    library                         => 4103684,
    mean_coverage                   => 66.88,
    mean_insert_size                => 314.4,
    paired                          => 1,
    platform                        => 'SLX',
    population                      => 'Population',
    project                         => 'SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants',
    reads                           => 11557928,
    reads_mapped                    => 11465632,
    reads_paired                    => 11434358,
    reverse_reads                   => 5778964,
    rmdup_bases                     => 1108486300,
    rmdup_bases_mapped              => 1097991197,
    rmdup_reads                     => 11084863,
    rmdup_reads_mapped              => 10992567,
    sample                          => 'SC_MFY5249244',
    sd_insert_size                  => 81.2,
    species                         => 'Pombe',
    study                           => 'ERP001017',
    targeted_avg_read_length        => 100,
    targeted_bases                  => 785559800,
    targeted_bases_mapped           => 783669700,
    targeted_bases_mapped_c         => 782746015,
    targeted_bases_of_100X_coverage => 236968,
    targeted_bases_of_10X_coverage  => 11374070,
    targeted_bases_of_1X_coverage   => 11418355,
    targeted_bases_of_20X_coverage  => 11207157,
    targeted_bases_of_2X_coverage   => 11417323,
    targeted_bases_of_50X_coverage  => 9395673,
    targeted_bases_of_5X_coverage   => 11409861,
    targeted_bases_trimmed          => 30414593,
    targeted_error_rate             => '3.641800e-03',
    targeted_filtered_reads         => 7855598,
    targeted_forward_reads          => 3927828,
    targeted_mean_coverage          => 66.13,
    targeted_mean_insert_size       => 316.3,
    targeted_paired                 => 1,
    targeted_reads                  => 7855598,
    targeted_reads_mapped           => 7836697,
    targeted_reads_paired           => 7817796,
    targeted_reverse_reads          => 3927770,
    targeted_rmdup_bases            => 777710400,
    targeted_rmdup_bases_mapped     => 774905867,
    targeted_rmdup_reads            => 7777104,
    targeted_rmdup_reads_mapped     => 7758203,
    targeted_sd_insert_size         => 80.9,
    withdrawn                       => 0
  },
  'metadata got updated for one of the improved bam files after redoing QC';

$passed_auto_qc_lanes = 0;
foreach my $bam (@improved_bams) {
    my $meta = VRPipe::File->get(path => $bam)->metadata;
    $passed_auto_qc_lanes++ if $meta->{auto_qc_status} eq 'passed';
}
is $passed_auto_qc_lanes, 5, 'after redoing QC, still have 5 lanes passed';

# mergeup
$output_dir = get_output_dir('lane_merge');
VRPipe::PipelineSetup->create(
    name       => 'pombe merge lanes',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'pombe improvement[10]',
        options => { metadata_keys => 'population|sample|platform|library' }
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'bam_merge_lanes_and_fix_rgs'),
    options     => {
        bam_tags_to_strip                     => 'OQ XM XG XO',
        bam_merge_keep_single_paired_separate => 1,
        cleanup                               => 1,
        delete_input_bams                     => 0
    }
);

$output_dir = get_output_dir('library_merge');
VRPipe::PipelineSetup->create(
    name       => 'pombe merge libraries',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'pombe merge lanes[4]',
        options => { metadata_keys => 'population|sample|platform' }
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'bam_merge_and_split'),
    options     => {
        bam_merge_keep_single_paired_separate => 0,
        split_bam_make_unmapped               => 0,
        delete_input_bams                     => 1,
        remove_merged_bams                    => 1
    }
);

# mergeacross
$output_dir = get_output_dir('merge_across');
my $mergeacross_ps = VRPipe::PipelineSetup->create(
    name       => 'pombe mergeacross',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'pombe merge libraries[2]',
        options => { metadata_keys => 'split_sequence' }
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'bam_merge'),
    options     => {
        bam_merge_keep_single_paired_separate => 0,
        bam_merge_maximum_files               => 2,
        delete_input_bams                     => 1
    }
);

handle_pipeline();
my @merged_bams;
foreach my $element (@{ get_elements($mergeacross_ps->datasource) }) {
    my @output_subdirs = output_subdirs($element->id, $mergeacross_ps->id);
    for my $chunk (1 .. 3) {
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
is_deeply \%seen_splits, { chromMT => 3, chromIII => 3, chromII => 3, chromI => 3, chromAB325691 => 3, chromMTR => 3 }, 'there were 3 mergeacross bams for each chromosome, and the metadata was correct';

finish;
