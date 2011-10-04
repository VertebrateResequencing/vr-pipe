#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;
use Cwd;

BEGIN {
    use Test::Most tests => 7;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $improvement_output_dir = get_output_dir('improvement_test');

ok my $improvement_pipeline = VRPipe::Pipeline->get(name => 'improvement_test_pipeline'), 'able to get the improvement_test_pipeline pipeline';

my @s_names;
foreach my $stepmember ($improvement_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(test_import_bams
                             sequence_dictionary
                             bam_index
                             gatk_target_interval_creator
                             bam_realignment_around_known_indels
                             bam_fix_mates
                             bam_index
                             bam_count_covariates
                             bam_recalibrate_quality_scores
                             bam_calculate_bq);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($improvement_output_dir, 'ref');
$improvement_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $res_dir = dir($improvement_output_dir, 'resources');
$improvement_pipeline->make_path($res_dir);

my $known_indels_source = file(qw(t data known_indels.vcf.gz));
my $known_indels = file($res_dir, 'known_indels.vcf.gz')->stringify;
copy($known_indels_source, $known_indels);
copy($known_indels_source.'.tbi', $known_indels.'.tbi');

my $known_sites_source = file(qw(t data known_sites.vcf.gz));
my $known_sites = file($res_dir, 'known_sites.vcf.gz')->stringify;
copy($known_sites_source, $known_sites);
copy($known_sites_source.'.tbi', $known_sites.'.tbi');

ok my $ds = VRPipe::DataSource->get(type => 'fofn',
                                 method => 'all',
                                 source => file(qw(t data improvement_datasource.fofn))->absolute->stringify,
                                 options => {}), 'could create a fofn datasource';

my @results = ();
while (my $element = $ds->next_element) {
    push(@results, $element->result);
}
my $cwd = cwd();
is_deeply \@results, [{paths => [file($cwd, 't', 'data', '2822_7.pe.bam')]}, 
                      {paths => [file($cwd, 't', 'data', '2822_6.pe.bam')]}, 
                      {paths => [file($cwd, 't', 'data', '2822_6.se.bam')]}, 
                      {paths => [file($cwd, 't', 'data', '2823_4.pe.bam')]}, 
                      {paths => [file($cwd, 't', 'data', '8324_8.pe.bam')]}], 'got correct results for fofn all';

my $improvement_pipelinesetup = VRPipe::PipelineSetup->get(name => 's_suis improvement',
                                                       datasource => VRPipe::DataSource->get(type => 'fofn',
                                                                                             method => 'all',
                                                                                             source => file(qw(t data improvement_datasource.fofn))->absolute->stringify,
                                                                                             options => {} ),
                                                       output_root => $improvement_output_dir,
                                                       pipeline => $improvement_pipeline,
                                                       options => {reference_fasta => $ref_fa,
                                                                   reference_assembly_name => 'SSuis1',
                                                                   reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                   reference_species => 'S.Suis',
                                                                   known_indels_for_realignment => "-known $known_indels",
                                                                   known_sites_for_recalibration => "-knownSites $known_sites",
                                                                   gatk_count_covariates_options => '-l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
                                                                   gatk_path => '/software/vertres/bin-external/GenomeAnalysisTK-1.2-29/',
                                                                   picard_path => '/software/vertres/bin-external/picard-tools-1.53/',
                                                                   cleanup => 0,
                                                                   sequence_dictionary_memory => 150,
                                                                   sequence_dictionary_time => 1});


my @output_files = (file($improvement_output_dir, 'ref', 'S_suis_P17.fa.dict'), file($improvement_output_dir, 'resources', 'known_indels.intervals'));

my @files = ('2822_7.pe.bam', '2822_6.pe.bam', '2822_6.se.bam', '2823_4.pe.bam', '8324_8.pe.bam');
my $element_id = 0;
foreach my $file (@files) {
    $element_id++;
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'test_import_bams', $file);
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'test_import_bams', "$file.bai");
    $file =~ s/bam$/realign.bam/;
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'bam_realignment_around_known_indels', $file);
    $file =~ s/bam$/sort.bam/;
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'bam_fix_mates', $file);
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'bam_fix_mates', "$file.bai");
    $file =~ s/bam$/recal_data.csv/;
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'bam_count_covariates', $file);
    $file =~ s/recal\_data\.csv$/recal.bam/;
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'bam_recalibrate_quality_scores', $file);
    $file =~ s/bam$/calmd.bam/;
    push @output_files, file($improvement_output_dir, output_subdirs($element_id), 'bam_calculate_bq', $file);
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

is_deeply [VRPipe::StepState->get(pipelinesetup => 1, stepmember => 4, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 5, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 6, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 8, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 9, dataelement => 1)->cmd_summary->summary,
           VRPipe::StepState->get(pipelinesetup => 1, stepmember => 10, dataelement => 1)->cmd_summary->summary],
          ['java $jvm_args -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference_fasta -o $intervals_file -known $known_indels_file(s) ',
           'java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing',
           'java $jvm_args -jar FixMateInformation.jar INPUT=$bam_file OUTPUT=$fixmate_bam_file SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0',
           'java $jvm_args -jar GenomeAnalysisTK.jar -T CountCovariates -R $reference_fasta -I $bam_file -recalFile $bam_file.recal_data.csv -knownSites $known_sites_file(s) -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
           'java $jvm_args -jar GenomeAnalysisTK.jar -T TableRecalibration -R $reference_fasta -recalFile $bam_file.recal_data.csv -I $bam_file -o $recalibrated_bam_file -l INFO --disable_bam_indexing',
           'samtools calmd -Erb $bam_file $reference_fasta > $bq_bam_file'],
          'cmd summaries for the major steps were as expected';
          
finish;