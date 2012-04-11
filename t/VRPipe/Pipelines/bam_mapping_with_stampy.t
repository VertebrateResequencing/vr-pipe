#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES',
                    required_exe => [qw(bwa stampy.py samtools)]);
    use TestPipelines;
}

my $mapping_output_dir = get_output_dir('bam_mapping_with_stampy');

ok my $mapping_pipeline = VRPipe::Pipeline->get(name => 'bam_mapping_with_stampy'), 'able to get the bam_mapping_with_stampy pipeline';

my @s_names;
foreach my $stepmember ($mapping_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(sequence_dictionary bwa_index stampy_buildgenome stampy_buildhash bam_metadata bam_name_sort bam_to_fastq fastq_split stampy_map_fastq sam_to_fixed_bam bam_merge_lane_splits)], 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($mapping_output_dir, 'ref');
$mapping_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);
my $mapping_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis mapping with stampy',
                                                       datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                             method => 'all',
                                                                                             source => file(qw(t data datasource.bam_fofn))),
                                                       output_root => $mapping_output_dir,
                                                       pipeline => $mapping_pipeline,
                                                       options => {fastq_chunk_size => 8000,
                                                                   reference_fasta => $ref_fa,
                                                                   reference_assembly_name => 'SSuis1',
                                                                   reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                   reference_species => 'S.Suis',
                                                                   bwa_index_options => '-a is',
                                                                   stampy_map_options => '--substitutionrate=0.002',
                                                                   stampy_bwa_options => '-q 15',
                                                                   cleanup => 1,
                                                                   sequence_dictionary_memory => 150,
                                                                   sequence_dictionary_time => 1,
                                                                   bwa_index_memory => 150,
                                                                   bwa_index_time => 1,
                                                                   bam_metadata_memory => 150,
                                                                   bam_metadata_time => 1,
                                                                   fastq_split_memory => 150,
                                                                   fastq_split_time => 1,
                                                                   sam_to_fixed_bam_memory => 150,
                                                                   sam_to_fixed_bam_time => 1,
                                                                   bam_merge_lane_splits_memory => 150,
                                                                   bam_merge_lane_splits_time => 1});

my $mapping_output_dir2 = get_output_dir('bam_mapping_with_stampy_divergent');

ok my $mapping_pipeline2 = VRPipe::Pipeline->get(name => 'bam_mapping_with_stampy_divergent'), 'able to get the bam_mapping_with_stampy_divergent pipeline';

@s_names = ();
foreach my $stepmember ($mapping_pipeline2->steps) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(sequence_dictionary stampy_buildgenome stampy_buildhash bam_metadata bam_name_sort bam_to_fastq fastq_split stampy_map_fastq sam_to_fixed_bam bam_merge_lane_splits bam_substitution_rate stampy_map_fastq sam_to_fixed_bam bam_merge_lane_splits bamcheck)], 'the pipeline has the correct steps';

VRPipe::PipelineSetup->get(name => 's_suis mapping with stampy',
                           datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                 method => 'all',
                                                                 source => file(qw(t data datasource.bam_fofn))),
                           output_root => $mapping_output_dir2,
                           pipeline => $mapping_pipeline2,
                           options => {fastq_chunk_size => 8000,
                                       reference_fasta => $ref_fa,
                                       reference_assembly_name => 'SSuis1',
                                       reference_public_url => 'ftp://s.suis.com/ref.fa',
                                       reference_species => 'S.Suis',
                                       stampy_map_options => '--substitutionrate=0.002',
                                       cleanup => 1,
                                       sequence_dictionary_memory => 150,
                                       sequence_dictionary_time => 1,
                                       bwa_index_memory => 150,
                                       bwa_index_time => 1,
                                       bam_metadata_memory => 150,
                                       bam_metadata_time => 1,
                                       fastq_split_memory => 150,
                                       fastq_split_time => 1,
                                       sam_to_fixed_bam_memory => 150,
                                       sam_to_fixed_bam_time => 1,
                                       bam_merge_lane_splits_memory => 150,
                                       bam_merge_lane_splits_time => 1});

my @bams;
my $element_num = 0;
foreach my $bam (qw(2822_6.se 2822_6.pe 2822_7.pe 2823_4.pe 8324_8.pe)) {
    $element_num++;
    my @output_subdirs = output_subdirs($element_num, 1);
    push(@bams, file(@output_subdirs, '11_bam_merge_lane_splits', "$bam.bam"));
    @output_subdirs = output_subdirs($element_num, 2);
    push(@bams, file(@output_subdirs, '14_bam_merge_lane_splits', "$bam.bam"));
}

ok handle_pipeline(@bams), 'pipelines ran ok and produced all bams';

my %scsss;
my $schema = $mapping_pipeline2->result_source->schema;
my $rs = $schema->resultset('StepCmdSummary');
while (my $scs = $rs->next) {
    my $summary = $scs->summary;
    next unless $summary =~ /^stampy\.py --substitutionrate/;
    $scsss{$summary} = 1;
}
is_deeply \%scsss, {'stampy.py --substitutionrate=0.002 --bwa=bwa --bwaoptions={-q 15 $ref} -g $ref.fa -h $ref.fa -o $out.sam -M $fastq(s)' => 1,
                    'stampy.py --substitutionrate=0.002 -g $ref.fa -h $ref.fa -o $out.sam -M $fastq(s)' => 1,
                    'stampy.py --substitutionrate=0.00001 -g $ref.fa -h $ref.fa -o $out.sam -M $fastq(s)' => 1}, 'the step command summarys for the stampy step were correct';

finish;