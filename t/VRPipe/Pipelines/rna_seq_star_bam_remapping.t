#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 10;
    use VRPipeTest;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER)],
        required_exe => [qw(imeta)]
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
}

# these are author-only tests (Sanger-specific) to make sure that we get all
# the desired metadata out of irods and warehouse for all the different types
# of data we process. We'll also test that we can then take this VRPipe metadata
# on files and store all the tracking information in an external database

my $output_root = get_output_dir('irods_metadata_populator');

# bam sequencing data
ok my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'seq',
    options => {
        file_query     => q[study_id = 2791 and sample = 'QC2Hip-1070' and target = 1],
        local_root_dir => $output_root
    }
  ),
  'could create an irods datasource for bam files';

my $pipeline = VRPipe::Pipeline->create(name => 'vrtrack_populate_from_irods_and_download_files');

create_fresh_vrtrack_test_db();
VRPipe::PipelineSetup->create(
    name        => 'bam populate',
    datasource  => $ds,
    output_root => $output_root,
    pipeline    => $pipeline,
    options     => {
        vrtrack_db                 => $ENV{VRPIPE_VRTRACK_TESTDB},
        vrlane_storage_dir         => $output_root,
        irods_download_input_files => 1
    }
);

my @expected_output_files = (file($output_root, '/seq/12820/12820_7#1.bam'), file($output_root, '/seq/12880/12880_6#1.bam'));
ok handle_pipeline(@expected_output_files), 'vrtrack_populate_from_vrpipe_metadata pipeline ran ok for bam files';

my $files = 0;
foreach my $result (map { result_with_inflated_paths($_) } @{ get_elements($ds) }) {
    $files++;
    next if $files > 1;
    my ($path) = @{ $result->{paths} };
    my $file = VRPipe::File->create(path => $path);
    is_deeply $file->metadata, {
        'alignment'               => '1',
        'id_run'                  => '12820',
        'irods_path'              => '/seq/12820/12820_7#1.bam',
        'is_paired_read'          => '1',
        'lane'                    => '12820_7#1',
        'library'                 => '10065808',
        'library_id'              => '10065808',
        'manual_qc'               => '1',
        'md5'                     => '7a5884b6b9f11c991ff856e27af7daf3',
        'public_name'             => 'HPSI0813i-fpdr',
        'reference'               => '/lustre/scratch110/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa/hs37d5.fa',
        'sample'                  => 'QC2Hip-1070',
        'sample_accession_number' => 'SAMEA2399188',
        'sample_cohort'           => '6d3d2acf-29a5-41a2-8992-1414706a527d',
        'sample_common_name'      => 'Homo sapiens',
        'sample_control'          => '0',
        'sample_created_date'     => '2014-02-17 13:40:22',
        'sample_id'               => '1860073',
        'sample_public_name'      => 'HPSI0813i-fpdr',
        'sample_supplier_name'    => 'face6d88-7e90-4215-aa80-fb2c3df5a4ed',
        'study'                   => 'G0325 [rnaseq] HipSci_RNASEQ_Normals_MA',
        'study_accession_number'  => 'EGAS00001000593',
        'study_id'                => '2791',
        'study_title'             => 'G0325 [rnaseq] HipSci_RNASEQ_Normals_MA',
        'tag'                     => 'ATCACG',
        'tag_index'               => '1',
        'target'                  => '1',
        'taxon_id'                => '9606',
        'total_reads'             => '24376496'
      
      },
      'correct file metadata was present on the first bam file';
}
is $files, 2, 'bam datasource returned the correct number of files';

my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
my @lanes = $vrtrack->get_lanes();
is @lanes, 2, 'the correct number of bam lanes were populated';
my $lane_info = lane_info(grep { $_->name eq '12820_7#1' } @lanes);

my $expected = {
    lane_name             => '12820_7#1',
    lane_hierarchy_name   => '12820_7#1',
    raw_reads             => 24376496,
    is_paired             => 1,
    lane_acc              => undef,
    storage_path          => undef,
    file_name             => '12820_7#1.bam',
    file_type             => 4,
    file_md5              => '7a5884b6b9f11c991ff856e27af7daf3',
    library_name          => 10065808,
    library_ssid          => 10065808,
    library_ts            => undef,
    project_id            => 2791,
    project_name          => 'G0325 [rnaseq] HipSci_RNASEQ_Normals_MA',
    study_acc             => 'EGAS00001000593',
    sample_name           => 'HPSI0813i-fpdr',
    sample_hierarchy_name => 'face6d88-7e90-4215-aa80-fb2c3df5a4ed',
    sample_ssid           => '1860073',
    species_taxon_id      => 9606,
    species_name          => 'Homo sapiens',
    individual_name       => '6d3d2acf-29a5-41a2-8992-1414706a527d',
    individual_acc        => 'SAMEA2399188'
};

is_deeply $lane_info, $expected, 'VRTrack was correctly populated for the first bam lane';

# create STAR setup using the output bam files from the import
my $output_dir = get_output_dir('rna_seq_star_bam_remapping');

# check pipeline has correct steps
ok my $star_pipeline = VRPipe::Pipeline->create(name => 'rna_seq_star_bam_remapping'), 'able to get the rna_seq_star_bam_remapping pipeline';
my @sb_names;
foreach my $stepmember ($star_pipeline->step_members) {
    push(@sb_names, $stepmember->step->name);
}

is_deeply \@sb_names, [qw(sequence_dictionary star_buildgenome bam_metadata bam_to_fastq star_map_fastq sam_to_fixed_bam bam_reheader bam_index)], 'the rna_seq_star_bam_remapping pipeline has the correct steps';

my $ref_fasta_orig = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($output_dir, 'ref');
$star_pipeline->make_path($ref_dir);
my $ref_fasta = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fasta_orig, $ref_fasta);

my $star_setup = VRPipe::PipelineSetup->create(
    name       => 'star setup',
    pipeline   => $star_pipeline,
    datasource => VRPipe::DataSource->create(
        type   => 'vrpipe',
        method => 'all',
        source => 'bam populate[2:input_files]',
    ),
    options => {
        reference_fasta => $ref_fasta,
    },
    output_root => $output_dir
);

my @final_files;
foreach my $name (qw(Genome SA SAindex)) {
    push(@final_files, file($ref_dir, $name));
}
foreach my $element (@{ get_elements($star_setup->datasource) }) {
    my @output_dirs = output_subdirs($element->id, $star_setup->id);
    push(@final_files, file(@output_dirs, '7_bam_reheader', 'Aligned.out.bam.bai'));
}
ok handle_pipeline(@final_files), 'rna_seq_star_bam_remapping pipeline ran ok';

finish;
exit;

sub create_fresh_vrtrack_test_db {
    my %args      = @_;
    my $add_notes = delete $args{with_notes};
    
    my %cd = VRTrack::Factory->connection_details('rw');
    open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
    print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    foreach my $sql (VRTrack::VRTrack->schema()) {
        print $mysqlfh $sql;
    }
    
    if ($add_notes) {
        print $mysqlfh q[insert into note set note_id = 1, note = "Control"; insert into note set note_id = 2, note = "Stem cell";];
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
        defined $control ? (control => $control) : ()
    };
    
    return $info;
}
