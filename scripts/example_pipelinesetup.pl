#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use VRPipe::Persistent::Schema;

my $pipeline = VRPipe::Pipeline->get(name => 'bam_mapping_with_bwa');

my $base_dir = '/lustre/scratch102/projects/mouse/Assembly/Remapping_Evaluation/';

my $datasource = VRPipe::DataSource->get(type => 'fofn',
                                         method => 'all',
                                         source => file($base_dir, 'nodbams.fofn')->stringify);

foreach my $ref_info (['sanger', 'NOD_Sanger_contigs_v1'],
                      ['oxford', 'NOD_Oxford_contigs_v1'], 
                      ['ebi', 'NOD_Ebi_contigs_v1']) {

    VRPipe::PipelineSetup->get(name => 'mouse bam remapping against the '.$ref_info->[1].' assembly',
                               datasource => $datasource,
                               output_root => dir($base_dir, 'mapping/'.$ref_info->[0].'/v1/NOD_Mouse_Genome/pipeline_output')->stringify,
                               pipeline => $pipeline,
                               options => { reference_fasta => file($base_dir, 'assemblies/'.$ref_info->[0].'/v1/NOD_Mouse_Genome.1kb.fa')->stringify,
                                            reference_assembly_name => $ref_info->[1],
                                            reference_species => 'Mouse',
                                            samtools_exe => '/software/vertres/bin-external/samtools-svn-981/samtools',
                                            bwa_exe => '/software/vertres/bin-external/bwa-0.5.9' });

}

exit;

