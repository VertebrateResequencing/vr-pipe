#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 4;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_filter_merge_and_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'vcf_filter_merge_and_annotate'), 'able to get the vcf_filter_merge_and_annotate pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(vcf_filter vcf_merge vcf_annotate vcf_consequences);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my vcf_filter_merge_and_annotate pipeline setup',
                                                    datasource => VRPipe::DataSource->get(type => 'delimited',
                                                                                          method => 'all_columns',
                                                                                          options => { delimiter => "\t" },
                                                                                          source => file(qw(t data datasource.vcfs))),
                                                    output_root => $output_dir,
                                                    pipeline => $pipeline,
                                                    options => {'vcf-annotate_options' => '-a /lustre/scratch106/projects/uk10k/ref/snps/g1k_dbsnp132_annot.tab.gz -d /lustre/scratch106/projects/uk10k/ref/snps/g1k_dbsnp132_annot_desc.txt -c CHROM,FROM,REF,ALT,-,-,INFO/KGPilot123,INFO/dbSNP ',
                                                                'vcf-filter_exe' => '/software/vertres/codebase/scripts/vcf-filter',
                                                                'vcf-filter_options' => '-f /nfs/vertres01/conf/uk10k_gatk_20110715.filter',
                                                                'vcf2consequences_options' => '-s Homo_sapiens -gerp -grantham',
                                                                cleanup => 0});

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('test1', 'test2') {
	$element_id++;
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vcf_merge', "${in}.filt.${in}a.filt.merged.vcf.gz"));
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vcf_merge', "${in}.filt.${in}a.filt.merged.vcf.gz.tbi"));
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vcf_annotate', "${in}.filt.${in}a.filt.merged.annot.vcf.gz"));
    push(@final_files, file($output_dir, output_subdirs($element_id), 'vcf_consequences', "${in}.filt.${in}a.filt.merged.annot.conseq.vcf.gz"));
    push(@final_files, file($output_dir, output_subdirs($element_id), 'vcf_consequences', "${in}.filt.${in}a.filt.merged.annot.conseq.vcf.gz.tbi"));
}
ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
