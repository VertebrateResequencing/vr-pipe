#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 21;
    use VRPipeTest;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER SAMTOOLS HTSLIB)],
        required_exe => [qw(imeta)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
    use_ok('VRPipe::Schema');
}

# these are author-only tests (Sanger-specific) to make sure that we can do qc
# based on files stored in irods

# we test sequencing_qc_from_irods_and_update_vrtrack pipeline here, which
# is the same as testing sequencing_qc_from_irods as well, since that just has
# one less step

my $output_root = get_output_dir('irods_qc');

ok my $pipeline = VRPipe::Pipeline->create(name => 'sequencing_qc_from_irods_and_update_vrtrack'), 'able to get the sequencing_qc_from_irods_and_update_vrtrack pipeline';
my @s_names;
foreach my $stepmember ($pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(samtools_fasta_gc_stats npg_cram_stats_parser plot_bamstats vrtrack_populate_from_graph_db);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

ok my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'seq',
    options => {
        file_query       => q[id_run = 16131 and target = 1 and type = cram],
        require_qc_files => 1,
    }
  ),
  'could create an irods datasource using all_with_warehouse_metadata method with require_qc_files option';

# artificially change the library on one of the lanes so that we can later test
# that the cram header is "wrong"
$ds->elements;
my $schema = VRPipe::Schema->create('VRTrack');
my $wrong_lane = $schema->get('Lane', { unique => '16131_7' });
my ($correct_lib) = $wrong_lane->related(incoming => { type => 'sequenced' });
my ($correct_sample) = $correct_lib->related(incoming => { type => 'prepared' });
my $wrong_lib = $schema->add('Library', { name => 'wrong_library_name', id => '999' });
$wrong_lib->relate_to($wrong_lane, 'sequenced', selfish => 1);
$correct_sample->relate_to($wrong_lib, 'prepared');

my $ref_fa_source = file(qw(t data GRCm38.MT.ref.fa));
my $ref_dir = dir($output_root, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

create_fresh_vrtrack_test_db();

VRPipe::PipelineSetup->create(
    name        => 'qc from irods',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db      => $ENV{VRPIPE_VRTRACK_TESTDB},
        reference_fasta => $ref_fa,
    }
);

my @expected_output_files;
my $element_id = 1;
foreach my $lane (1 .. 8) {
    my @output_subdirs = output_subdirs($element_id++);
    foreach my $kind (qw(quals-hm quals quals2 quals3 insert-size gc-content gc-depth acgt-cycles coverage indel-dist indel-cycles)) {
        push(@expected_output_files, file(@output_subdirs, '3_plot_bamstats', "16131_$lane" . '-' . $kind . '.png'));
    }
}

ok handle_pipeline(@expected_output_files), 'sequencing_qc_from_irods pipeline downloaded all the expected qc files';

my $files = 0;
foreach my $result (map { result_with_inflated_files($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($file) = @{ $result->{files} };
    $file = $schema->get_file($file->protocolless_path, $file->protocol);
    my $props = $file->properties(flatten_parents => 1);
    delete @{$props}{qw(datasource_uuid uuid filesystemelement_uuid stepstate_uuid dataelement_uuid stepstate_sql_id datasource_sql_id dataelement_sql_id)};
    is_deeply $props,
      {
        'vrtrack_lane_run'             => '16131',
        'vrtrack_library_id'           => '13607695',
        'vrtrack_sample_supplier_name' => 'XX339024',
        'vrtrack_donor_id'             => 'vbseqx106041566',
        'target'                       => '1',
        'vrtrack_sample_control'       => '0',
        'vrtrack_alignment_reference'  => '/lustre/scratch110/srpipe/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/bwa0_6/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa',
        'datasource_method'            => 'all_with_warehouse_metadata',
        'vrtrack_study_id'             => '3474',
        'vrtrack_taxon_common_name'    => 'Homo Sapien',
        'vrtrack_sample_name'          => 'vbseqx106041566',
        'vrtrack_sample_id'            => '2218372',
        'manual_qc'                    => '1',
        'vrtrack_sample_created_date'  => '1422355113',
        'vrtrack_lane_total_reads'     => '747630612',
        'vrtrack_library_name'         => '13607695',
        'filesystemelement_path'       => '/seq',
        'md5'                          => '67de7687a3c87570b0c28af67bf3f7f5',
        'vrtrack_study_name'           => 'IHTP_MWGS ESGI - WGS of samples from the INGI-Val Borbera genetic isolate (X10)',
        'vrtrack_lane_is_paired_read'  => '1',
        'vrtrack_taxon_id'             => '9606',
        'vrtrack_study_accession'      => 'EGAS00001001123',
        'vrtrack_sample_accession'     => 'EGAN00001264662',
        'vrtrack_lane_lane'            => '1',
        'vrtrack_group_name'           => 'all_studies',
        'pipelinesetup_id'             => '1',
        'vrtrack_lane_unique'          => '16131_1',
        'filesystemelement_basename'   => 'irods:/',
        'path'                         => '/seq/16131/16131_1.cram',
        'datasource_source'            => 'seq',
        'pipelinesetup_name'           => 'qc from irods',
        'basename'                     => '16131_1.cram',
        'pipelinesetup_output_root'    => $output_root,
        'datasource_type'              => 'irods',
        'pipelinesetup_user'           => scalar(getpwuid($<))
      },
      'the graph database had the correct properties on nodes related to the first cram file';
}
is $files, 8, 'irods datasource returned the correct number of files';

ok my $lane = $schema->get('Lane', { unique => '16131_8' }), 'One of the lanes was in the db';
my ($stats)    = $lane->related(outgoing => { max_depth => 5, namespace => 'VRTrack', label => 'Bam_Stats' });
my ($geno)     = $lane->related(outgoing => { max_depth => 5, namespace => 'VRTrack', label => 'Genotype' });
my ($verify)   = $lane->related(outgoing => { max_depth => 5, namespace => 'VRTrack', label => 'Verify_Bam_ID' });
my ($mistakes) = $lane->related(outgoing => { max_depth => 5, namespace => 'VRTrack', label => 'Header_Mistakes' });
my $all_made = $stats && $geno && $verify && $mistakes;
ok $all_made, 'There were related Bam_Stats, Genotype, Verify_Bam_ID and Header_Mistakes nodes in the db';

my $props = $stats->properties;
delete @{$props}{qw(date uuid)};
is_deeply $props,
  {
    'bases trimmed'                  => '0',
    'average quality'                => '36.0',
    'mode'                           => 'normal',
    'is sorted'                      => '1',
    'reads mapped after rmdup'       => '650499184',
    'mismatches'                     => '1150869760',
    'reads duplicated'               => '70304307',
    'bases of 2X coverage'           => '2909978551',
    'bases of 10X coverage'          => '2894225756',
    'bases of 5X coverage'           => '2903606279',
    'raw total sequences'            => '813763252',
    'reads paired'                   => '721009154',
    'maximum length'                 => '151',
    'insert size average'            => '531.3',
    'non-primary alignments'         => '0',
    'reads mapped and paired'        => '720613244',
    'reads after rmdup'              => '743458945',
    'bases of 1X coverage'           => '2914630840',
    'bases duplicated'               => '10615950357',
    'bases of 50X coverage'          => '142529436',
    'average length'                 => '151',
    '1st fragments'                  => '360504577',
    'reads properly paired'          => '690246846',
    'options'                        => '-F 0xB00',
    'outward oriented pairs'         => '2803983',
    'bases of 100X coverage'         => '13659010',
    'bases mapped (cigar)'           => '105658938272',
    'reads unmapped'                 => '205663',
    'bases mapped'                   => '108841327141',
    'filtered sequences'             => '92754098',
    'inward oriented pairs'          => '346653690',
    'pairs with other orientation'   => '950678',
    'sequences'                      => '721009154',
    'bases after rmdup'              => '98256431897',
    'reads MQ0'                      => '42278001',
    'bases mapped after rmdup'       => '95353385474',
    'insert size standard deviation' => '782.8',
    'pairs on different chromosomes' => '9355186',
    'mean coverage'                  => '36.18',
    'reads mapped'                   => '720803491',
    'bases of 20X coverage'          => '2831811694',
    'error rate'                     => '1.089231e-02',
    'last fragments'                 => '360504577',
    'reads QC failed'                => '0',
    'total length'                   => '108872382254'
  },
  'Bam_Stats node had the correct properties';

$props = $geno->properties;
delete @{$props}{qw(date uuid)};
is_deeply $props,
  {
    expected_sample_name     => 'vbseqx106041569',
    pass                     => 1,
    bam_gt_likelihood_string => '122,0,255;0;237,99,0;0;203,60,0;255,0,249;0;201,63,0;243,75,0;124,0,249;0;0;255,66,0;178,0,255;202,81,0;0;255,44,0;0;241,105,0;0;255,0,255;255,66,0;84,0,255;233,57,0;248,0,255;191,0,223',
    bam_gt_depths_string     => '33;34;35;29;23;29;30;26;30;32;30;41;25;29;30;29;25;29;36;34;30;26;30;20;37;22',
    bam_call_count           => 26,
    bam_call_string          => 'TGAATTAACCCTAACCAATGTTGGAAAGAACCTTCCTTGGTCGGGACCGATC',
    matched_sample_name      => 'vbseqx106041569',
    match_score              => '1;0;1;1;1;1;1;0;1;1;1;1;1;1;1;0;1;0;1;1;1;1;1;1;1;1',
    common_snp_count         => 22,
    match_count              => 22,
    mismatch_count           => 0
  },
  'Genotype node had the correct properties';

$props = $verify->properties;
delete $props->{date};
delete $props->{uuid};
is_deeply $props,
  {
    freemix        => '0.00087',
    pass           => 1,
    number_of_snps => 1347021,
    avg_depth      => 29.79,
    freeLK0        => 9404857.63,
    freeLK1        => 9399169.37
  },
  'Verify_Bam_ID node had the correct properties';

$props = $mistakes->properties;
delete $props->{uuid};
is_deeply $props, { num_mistakes => 0, md5_of_ref_seq_md5s => 'f95dcc1c1300f59b028fba79f49878a4' }, 'Header_Mistakes node had the correct properties';

# and test the one we artificially made a mistake for
my ($fake_mistakes) = $wrong_lane->related(outgoing => { max_depth => 5, namespace => 'VRTrack', label => 'Header_Mistakes' });
$props = $fake_mistakes->properties;
delete $props->{uuid};
is_deeply $props,
  {
    num_mistakes        => 1,
    md5_of_ref_seq_md5s => 'f95dcc1c1300f59b028fba79f49878a4',
    LB                  => ['13607731', '999']
  },
  'Header_Mistakes nodes can correctly show mistakes';

my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
my @lanes = $vrtrack->get_lanes();
is @lanes, 8, 'the correct number of lanes were populated in the vrtrack db';
my ($first_lane) = grep { $_->name eq '16131_1' } @lanes;
my $lane_info = lane_info($first_lane);

my $expected = {
    'species_name'          => 'Homo Sapien',
    'library_ssid'          => '13607695',
    'library_ts'            => undef,
    'is_paired'             => 1,
    'file_name'             => '16131_1.cram',
    'individual_alias'      => '',
    'sample_name'           => 'vbseqx106041566',
    'project_name'          => 'IHTP_MWGS ESGI - WGS of samples from the INGI-Val Borbera genetic isolate (X10)',
    'lane_name'             => '16131_1',
    'project_id'            => '3474',
    'file_md5'              => '67de7687a3c87570b0c28af67bf3f7f5',
    'individual_name'       => 'vbseqx106041566',
    'sample_ssid'           => '2218372',
    'lane_acc'              => undef,
    'sample_hierarchy_name' => 'XX339024',
    'species_taxon_id'      => '9606',
    'library_name'          => '13607695',
    'study_acc'             => 'EGAS00001001123',
    'lane_hierarchy_name'   => '16131_1',
    'file_type'             => '6',
    'npg_qc_status'         => 'pass',
    'individual_acc'        => 'EGAN00001264662',
    'raw_reads'             => '724197352',
    'storage_path'          => undef
};

is_deeply $lane_info, $expected, 'VRTrack was correctly populated for the first cram lane';

ok my $mapstats = $first_lane->latest_mapping, 'the first lane had a mapstats in the VRTrack db';
is_deeply [$mapstats->raw_reads, $mapstats->clip_bases, $mapstats->reads_mapped, $mapstats->error_rate], [724197352, 109353800152, 723945472, 0.0113805], 'the mapstats had some correct values stored'; #*** assuming 'clip_bases' is supposed to return number of bases after clipping, and not the number of clipped bases
ok my $images = $mapstats->images, 'the mapstats had images';
is_deeply [sort map { $_->caption } @$images], ['ACGT Cycles', 'Coverage', 'GC Content', 'GC Depth', 'Indel distribution', 'Indels per cycle', 'Insert Size', 'Qualities', 'Qualities', 'Qualities', 'Qualities'], 'the images had the correct captions';

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

sub lane_info {
    my $lane   = shift;
    my ($file) = @{ $lane->files };
    my %h      = $vrtrack->lane_hierarchy_objects($lane);
    
    my $control;
    my $note_id = $h{sample}->note_id;
    if (defined $note_id) {
        $control = $note_id == 1 ? 1 : 0;
    }
    
    my $info = {
        lane_name             => $lane->name,
        lane_hierarchy_name   => $lane->hierarchy_name,
        raw_reads             => $lane->raw_reads,
        is_paired             => $lane->is_paired,
        lane_acc              => $lane->acc,
        storage_path          => $lane->storage_path,
        npg_qc_status         => $lane->npg_qc_status,
        file_name             => $file->name,
        file_type             => $file->type,
        file_md5              => $file->md5,
        library_name          => $h{library}->name,
        library_ssid          => $h{library}->ssid,
        library_ts            => $h{library}->library_tag_sequence,
        project_id            => $h{project}->ssid,
        project_name          => $h{project}->name,
        study_acc             => $h{study}->acc,
        sample_name           => $h{sample}->name,
        sample_hierarchy_name => $h{sample}->hierarchy_name,
        sample_ssid           => $h{sample}->ssid,
        species_taxon_id      => $h{species}->taxon_id,
        species_name          => $h{species}->name,
        individual_name       => $h{individual}->name,
        individual_acc        => $h{individual}->acc,
        individual_alias      => $h{individual}->alias,
        defined $control ? (control => $control) : ()
    };
    
    return $info;
}
