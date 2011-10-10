#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    
    use_ok('VRPipe::Persistent::Schema');
    
    use TestPipelines;
}

my $output_dir = get_output_dir('vcf_filter_merge_and_vep_annotate_pipeline');

ok my $pipeline = VRPipe::Pipeline->get(name => 'vcf_filter_merge_and_vep_annotate'), 'able to get the vcf_filter_merge_and_vep_annotate pipeline';
my @s_names;
foreach my $stepmember ($pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(vcf_filter vcf_merge vcf_annotate vep_analysis vcf_vep_consequences);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $test_pipelinesetup = VRPipe::PipelineSetup->get(name => 'my vcf_filter_merge_and_vep_annotate pipeline setup',
                                                    datasource => VRPipe::DataSource->get(type => 'delimited',
                                                                                          method => 'all_columns',
                                                                                          options => { delimiter => "\t" },
                                                                                          source => file(qw(t data datasource.vcfs))),
                                                    output_root => $output_dir,
                                                    pipeline => $pipeline,
                                                    options => {'vcf-annotate_options' => '-a /lustre/scratch106/projects/uk10k/ref/snps/g1k_dbsnp132_annot.tab.gz -d /lustre/scratch106/projects/uk10k/ref/snps/g1k_dbsnp132_annot_desc.txt -c CHROM,FROM,REF,ALT,-,-,INFO/KGPilot123,INFO/dbSNP ',
                                                                'vep_options' => '--sift b --polyphen b --condel b --gene --hgnc --format vcf --force_overwrite --cache --dir /lustre/scratch106/user/cj5/vep_cache',
                                                                'vep_exe' => '/nfs/users/nfs_c/cj5/vr-codebase/scripts/variant_effect_predictor.pl',
                                                                'vcf2consequences_exe' => '/nfs/users/nfs_c/cj5/vr-codebase/scripts/vcf2consequences_vep',
                                                                'vcf2consequences_options' => '-grantham --gerp /lustre/scratch106/user/cj5/gerp_db/hs',
                                                                cleanup => 0});

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('test1.vcf', 'test2.vcf') {
	$element_id++;
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vcf_merge', $in.'.filtered.merged.gz'));
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vcf_merge', $in.'.filtered.merged.gz.tbi'));
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vcf_annotate', $in.'.filtered.merged.annotated.gz'));
    push(@output_files, file($output_dir, output_subdirs($element_id), 'vep_analysis', $in.'.filtered.merged.annotated.vep.txt'));
    push(@final_files, file($output_dir, output_subdirs($element_id), 'vcf_consequences', $in.'.filtered.merged.annotated.conseq.gz'));
    push(@final_files, file($output_dir, output_subdirs($element_id), 'vcf_consequences', $in.'.filtered.merged.annotated.conseq.gz.tbi'));
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

done_testing;
exit;
