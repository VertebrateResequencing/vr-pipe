#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Path::Class qw(file dir);

BEGIN {
    use Test::Most tests => 5;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPersistentReal;
}

my $scheduler = VRPipe::Scheduler->get();
my $mapping_output_dir = dir($scheduler->output_root, 'pipelines_test_output', 'mapping');
$scheduler->remove_tree($mapping_output_dir);
$scheduler->make_path($mapping_output_dir);
my $manager = VRPipe::Manager->get();


ok my $mapping_pipeline = VRPipe::Pipeline->get(name => 'mapping'), 'able to get a pre-written pipeline';
TODO: {
    local $TODO = "not all steps implemented yet";
    my @s_names;
    foreach my $stepmember ($mapping_pipeline->steps) {
        push(@s_names, $stepmember->step->name);
    }
    is_deeply \@s_names, [qw(fastq_metadata fastq_split fastq_map bam_merge bam_stats)], 'the pipeline has the correct steps';
    
    $mapping_pipeline = VRPipe::Pipeline->get(name => 'mapping');
    @s_names = ();
    foreach my $stepmember ($mapping_pipeline->steps) {
        push(@s_names, $stepmember->step->name);
    }
    is_deeply \@s_names, [qw(fastq_metadata fastq_split fastq_map bam_merge bam_stats)], 'the pipeline has the correct steps after a second retrieval';
}

my $ref_fa = file(qw(t data S_suis_P17.fa));
my $mapping_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis mapping',
                                                       datasource => VRPipe::DataSource->get(type => 'sequence_index',
                                                                                             method => 'lane_fastqs',
                                                                                             source => file(qw(t data datasource.sequence_index))),
                                                       output_root => $mapping_output_dir,
                                                       pipeline => $mapping_pipeline);

my @mapping_output_files;

handle_pipeline(@mapping_output_files);

is_deeply [VRPipe::File->get(path => file(qw(t data 8324_8_1.fastq))->absolute)->md5,
           VRPipe::File->get(path => file(qw(t data 8324_8_1.fastq))->absolute)->metadata,
           VRPipe::File->get(path => file(qw(t data 8324_8_2.fastq))->absolute)->metadata],
          ['b90eb321fb12b47800132b237a6a3722',
           { lane => '8324_8',
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID02',
             sample => 'SAMPLE02',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB03',
             insert_size => 200,
             withdrawn => 0,
             reads => 250,
             bases => 15250,
             analysis_group => 'low coverage',
             paired => 1,
             mate => file(qw(t data 8324_8_2.fastq))->absolute->stringify },
           { expected_md5 => 'c8ecf1168b710adde2a3b45f98849780',
             lane => '8324_8',
             study => 'STUDY01',
             study_name => 'my study name',
             center_name => 'SC',
             sample_id => 'SAMPLEID02',
             sample => 'SAMPLE02',
             population => 'POP',
             platform => 'ILLUMINA',
             library => 'LIB03',
             insert_size => 200,
             withdrawn => 0,
             reads => 250,
             bases => 13500,
             analysis_group => 'low coverage',
             paired => 2,
             mate => file(qw(t data 8324_8_1.fastq))->absolute->stringify }], 'fastqs that went through the first step of the mapping pipeline have the correct metadata';

#is handle_pipeline(@mapping_output_files), 1, 'all mapping files were created via Manager';

done_testing;
exit;

sub handle_pipeline {
    my $give_up = 200;
    while (! $manager->trigger) {
        last if $give_up-- <= 0;
        $manager->handle_submissions;
        sleep(1);
    }
    my $all_created = 1;
    foreach my $ofile (@_) {
        unless (-s $ofile) {
            warn "$ofile is missing\n";
            $all_created = 0;
        }
    }
    return $all_created;
}