#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES CRAMTOOLS)],
                    required_exe => [qw(samtools bamcheck)]);
    use TestPipelines;
}

my $output_dir = get_output_dir('cram_test');

ok my $pipeline = VRPipe::Pipeline->create(name => 'bam_cram_bam_test_pipeline'), 'able to get the bam_cram_bam_test_pipeline pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(test_import_bams fasta_index bam_metadata bam_index bam_to_cram cram_index cram_to_bam)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $pipelinesetup = VRPipe::PipelineSetup->create(name       => 'bam_to_cram_to_bam',
                                                  datasource => VRPipe::DataSource->create(type    => 'fofn',
                                                                                           method  => 'all',
                                                                                           source  => file(qw(t data hs_chr20.qc.bam.fofn)),
                                                                                           options => {}),
                                                  output_root => $output_dir,
                                                  pipeline    => $pipeline,
                                                  options     => { reference_fasta => $ref_fa });

my @output_files = (file($output_dir, 'ref', 'human_g1k_v37.chr20.fa.fai'));

my @files = ('hs_chr20.a.bam', 'hs_chr20.b.bam', 'hs_chr20.c.bam', 'hs_chr20.d.bam');
my $element_id = 0;
foreach my $file (@files) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push @output_files, file(@output_subdirs, '1_test_import_bams', $file);
    push @output_files, file(@output_subdirs, '1_test_import_bams', "$file.bai");
    $file =~ s/bam$/cram/;
    push @output_files, file(@output_subdirs, '5_bam_to_cram', $file);
    push @output_files, file(@output_subdirs, '5_bam_to_cram', "$file.crai");
    $file =~ s/cram$/bam/;
    push @output_files, file(@output_subdirs, '7_cram_to_bam', $file);
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

finish;
