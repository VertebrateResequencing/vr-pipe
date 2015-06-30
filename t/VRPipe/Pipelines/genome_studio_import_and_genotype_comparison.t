#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;

BEGIN {
    use Test::Most tests => 40;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_AUTHOR_TESTS WAREHOUSE_DATABASE WAREHOUSE_HOST WAREHOUSE_PORT WAREHOUSE_USER)],
        required_exe => [qw(iget iquest fcr-to-vcf sort bgzip bcftools)]
    );
    use TestPipelines;
    use_ok('VRPipe::Schema');
}

my $output_dir = get_output_dir('genome_studio');
my $irods_dir = dir($output_dir, 'irods_import')->stringify;

# create a setup to import idat files from irods
my $ds = VRPipe::DataSource->create(
    type    => 'irods',
    method  => 'all_with_warehouse_metadata',
    source  => 'archive',
    options => {
        file_query     => q[study_id = 2624 and type = gtc and beadchip = 9439653037],
        local_root_dir => $irods_dir
    }
);

my $pipeline = VRPipe::Pipeline->create(name => 'irods_analysis_files_download');

my $setup1 = VRPipe::PipelineSetup->create(
    name        => 'gtc import',
    datasource  => $ds,
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        vrlane_storage_dir         => $irods_dir,
        irods_download_input_files => 1
    }
);

my @analysis_files = (file($irods_dir, '/archive/GAPI/gen/analysis/de/a3/63/coreex_hips/20140620/coreex_hips_20140620.fcr.txt.gz'));
ok handle_pipeline(@analysis_files), 'irods_analysis_files_download pipeline ran ok and got the analysis files';
$setup1->active(0);
$setup1->update;

# create split->vcf setup using the output gtc files from the import
$output_dir = get_output_dir('genome_studio_split_and_convert_to_vcf');

# check pipeline has correct steps
ok my $split_convert_pipeline = VRPipe::Pipeline->create(name => 'genome_studio_split_and_convert_to_vcf'), 'able to get the genome_studio_split_and_convert_to_vcf pipeline';
my @s_names;
foreach my $stepmember ($split_convert_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(split_genome_studio_genotype_files illumina_coreexome_manifest_to_map genome_studio_fcr_to_vcf)], 'the pipeline has the correct steps';

my $external_reheader_penncnv = file(qw(t data reheader_penncnv.txt))->absolute->stringify;
my $manifest_file_soure       = file(qw(t data hipsci_genotyping.snp.manifest))->absolute->stringify;
my $manifest_file             = file($output_dir, 'hipsci_genotyping.manifest')->stringify;
copy($manifest_file_soure, $manifest_file);

# create pipeline setup
my $split_convert_setup = VRPipe::PipelineSetup->create(
    name       => 'gtc split and convert',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'gtc import[1:input_files]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $split_convert_pipeline,
    options     => {
        reheader_penncnv   => $external_reheader_penncnv,
        coreexome_manifest => $manifest_file,
        cleanup            => 0
    }
);

my @samples    = qw(HPSI1013i-garx_2_QC1Hip-2193 HPSI1013pf-garx_QC1Hip-2191 HPSI1013i-funy_1_QC1Hip-2196 HPSI1013i-garx_1_QC1Hip-2192 HPSI0813pf-uaqe_QC1Hip-2188 HPSI1013pf-funy_QC1Hip-2195 HPSI0813i-uaqe_2_QC1Hip-2190 HPSI0613i-dium_3_QC1Hip-2224 HPSI1013i-funy_2_QC1Hip-2197 HPSI0813i-uaqe_1_QC1Hip-2189 HPSI1013i-garx_3_QC1Hip-2194 HPSI1013i-funy_3_QC1Hip-2198);
my $element_id = 12;
my @genotype_files;
my @vcf_files;
my @vcf_files_with_control;
foreach my $sample (@samples) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 2);
    
    # for testing purposes we'll fake a sample swap between 2 of the cohorts by
    # altering sample names and cohort
    if ($sample eq 'HPSI1013i-garx_1_QC1Hip-2192' || $sample eq 'HPSI1013i-funy_2_QC1Hip-2197') {
        my ($input_file) = @{ VRPipe::DataElement->get(id => $element_id)->files };
        if ($sample eq 'HPSI1013i-garx_1_QC1Hip-2192') {
            $input_file->add_metadata({ public_name => 'HPSI1013i-funy_2', sample => 'QC1Hip-2197', sample_cohort => '446ed213-1c31-427b-b314-2636ad8188d5' });
            $sample = 'HPSI1013i-funy_2_QC1Hip-2197';
        }
        else {
            $input_file->add_metadata({ public_name => 'HPSI1013i-garx_1', sample => 'QC1Hip-2192', sample_cohort => 'a34ff157-ce3c-46a1-b1e3-9349cfe8cd86' });
            $sample = 'HPSI1013i-garx_1_QC1Hip-2192';
        }
    }
    
    my ($sanger_name) = $sample =~ /(QC\S+)/;
    push(@genotype_files, file(@output_subdirs, '1_split_genome_studio_genotype_files', $sanger_name . '.genotyping.fcr.txt'));
    
    my $vcf_file = file(@output_subdirs, '3_genome_studio_fcr_to_vcf', "$sample/$sample.vcf.gz");
    push(@vcf_files, $vcf_file);
    push(@vcf_files_with_control, $vcf_file) unless ($element_id == 20);
}

# run pipeline and check outputs
ok handle_pipeline(@genotype_files, @vcf_files), 'genome_studio_split_and_convert_to_vcf pipeline ran ok';

# check genotype file metadata
my $meta = VRPipe::File->get(path => $genotype_files[0])->metadata;
is_deeply $meta,
  {
    analysis_uuid           => 'c837da60-20e7-43a5-bb52-b24759a4033d',
    beadchip                => 9439653037,
    beadchip_section        => 'R06C01',
    expected_md5            => '00d2cd3e5df1fd5fc73476458a75107b',
    infinium_plate          => 'WG0207474-DNA',
    infinium_sample         => '344265_F04_QC1Hip-2193',
    infinium_well           => 'F04',
    irods_analysis_files    => [qw(/archive/GAPI/gen/analysis/de/a3/63/coreex_hips/20140620/coreex_hips_20140620.fcr.txt.gz  /archive/GAPI/gen/analysis/de/a3/63/coreex_hips/20140620/genotyping.db)],
    irods_local_storage_dir => $irods_dir,
    irods_path              => '/archive/GAPI/gen/infinium/00/d2/cd/9439653037_R06C01.gtc',
    md5                     => '00d2cd3e5df1fd5fc73476458a75107b',
    public_name             => 'HPSI1013i-garx_2',
    sample                  => 'QC1Hip-2193',
    sample_accession_number => 'SAMEA2398552',
    sample_cohort           => 'a34ff157-ce3c-46a1-b1e3-9349cfe8cd86',
    sample_common_name      => 'Homo sapiens',
    sample_consent          => 1,
    sample_control          => 0,
    sample_created_date     => '2014-05-08 15:33:27',
    sample_donor_id         => 'a34ff157-ce3c-46a1-b1e3-9349cfe8cd86',
    sample_id               => 1943009,
    sample_supplier_name    => 'b91c722d-ee69-46db-bd02-5a32867bb838',
    study_id                => 2624,
    study_title             => 'G0325 [coreex] Wellcome Trust Strategic Award application – HIPS',
    taxon_id                => 9606
  },
  'metadata correct for one of the genotype files';

$meta = VRPipe::File->get(path => $vcf_files[0])->metadata;
is $meta->{sample_cohort}, 'a34ff157-ce3c-46a1-b1e3-9349cfe8cd86', 'the VCF file has sample_cohort metadata';

# we'll take this opportunity to test the vcf_merge_and_compare_genotypes
# pipeline as well
$output_dir = get_output_dir('vcf_merge_and_compare_genotypes');
my $vrpipe_ds = VRPipe::DataSource->create(
    type    => 'vrpipe',
    method  => 'group_by_metadata',
    source  => '2[3:vcf_files]',
    options => { metadata_keys => 'sample_cohort' }
);

# check pipeline has correct steps
ok my $gt_pipeline = VRPipe::Pipeline->create(name => 'vcf_merge_and_compare_genotypes'), 'able to get the vcf_merge_and_compare_genotypes pipeline';
@s_names = ();
foreach my $stepmember ($gt_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
is_deeply \@s_names, [qw(vcf_merge_different_samples vcf_genotype_comparison)], 'the pipeline has the correct steps';

# create pipeline setup
VRPipe::PipelineSetup->create(
    name        => 'genotype comparison',
    datasource  => $vrpipe_ds,
    output_root => $output_dir,
    pipeline    => $gt_pipeline,
    options     => {}
);

my (@merged_vcf_files, @gtypex_files, @expected_metadata);
# check whether swapped samples can be identified by genotype comparision pipeline
foreach my $element_id (25 .. 28) {
    my @output_subdirs = output_subdirs($element_id, 3);
    push(@merged_vcf_files, file(@output_subdirs, '1_vcf_merge_different_samples', 'merged.vcf.gz'));
    push(@gtypex_files,     file(@output_subdirs, '2_vcf_genotype_comparison',     'merged.vcf.gz.gtypex'));
    
    my $group = VRPipe::DataElement->get(id => $element_id)->metadata->{group};
    my %expected = (group => $group);
    if ($group eq '9f6e6b87-957a-4226-b9b9-ea0734d1b8c1') {
        $expected{genotype_maximum_deviation} = ['==', 0, 'HPSI0613i-dium_3_QC1Hip-2224'];
    }
    elsif ($group eq '750c3684-c29d-4ead-b8fb-3d9ded319b83') {
        $expected{genotype_maximum_deviation} = ['>=', 1, 'HPSI0813i-uaqe_1_QC1Hip-2189'];
    }
    elsif ($group eq 'a34ff157-ce3c-46a1-b1e3-9349cfe8cd86') {
        $expected{genotype_maximum_deviation} = ['>=', 10, 'HPSI1013i-garx_1_QC1Hip-2192'];
    }
    elsif ($group eq '446ed213-1c31-427b-b314-2636ad8188d5') {
        $expected{genotype_maximum_deviation} = ['>=', 10, 'HPSI1013i-funy_2_QC1Hip-2197'];
    }
    push(@expected_metadata, \%expected);
}

ok handle_pipeline(@merged_vcf_files, @gtypex_files), 'vcf_merge_and_compare_genotypes pipeline ran ok';

foreach my $vcf_path (@merged_vcf_files) {
    $meta = VRPipe::File->get(path => $vcf_path)->metadata;
    my $expected = shift @expected_metadata;
    
    is "$meta->{sample_cohort}", $expected->{group}, "sample_cohort and beadchip metadata was correct for one of the merged VCF files";
    
    my ($cmp, $eval, $esample) = @{ $expected->{genotype_maximum_deviation} };
    my ($aval, $asample) = split(':', $meta->{genotype_maximum_deviation});
    cmp_ok $aval, $cmp, $eval, "genotype_maximum_deviation metadata value was correct for one of the merged VCF files";
    is $asample, $esample, "genotype_maximum_deviation metadata sample was correct for one of the merged VCF files";
}

# create fofn datasource for cnv and loh pipelines
my $fofn_file = file($output_dir, 'sample_vcfs.fofn');
my $fh = $fofn_file->openw;
print $fh "path\tsample_cohort\n";
foreach my $path (@vcf_files_with_control) {
    my $file = VRPipe::File->get(path => $path);
    my $sample_cohort = $file->metadata->{sample_cohort};
    print $fh "$path\t$sample_cohort\n";
}
close($fh);
my $fofn_ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'grouped_by_metadata',
    options => { metadata_keys => 'sample_cohort' },
    source  => $fofn_file->stringify
);

# we'll also take the opportunity to test hipsci cnv caller pipeline, since that
# also uses files from the genome studio import
SKIP: {
    my $num_tests = 10;
    skip "hipsci cnv calling tests disabled without plot-polysomy.py and cmp-cnvs.pl in your path", $num_tests unless can_execute('bcftools') && can_execute('cmp-cnvs.pl') && can_execute('plot-polysomy.py');
    
    $output_dir = get_output_dir('polysomy_cnv_caller');
    
    # check pipeline has correct steps
    ok my $cnv_pipeline = VRPipe::Pipeline->create(name => 'polysomy_cnv_caller'), 'able to get the polysomy_cnv_caller pipeline';
    my @step_names;
    foreach my $stepmember ($cnv_pipeline->step_members) {
        push(@step_names, $stepmember->step->name);
    }
    is_deeply \@step_names, [qw(vcf_merge_different_samples_control_aware polysomy plot_polysomy bcftools_cnv combine_bcftools_cnvs)], 'the polysomy_cnv_caller pipeline has the correct steps';
    
    my $cnv_setup = VRPipe::PipelineSetup->create(
        name        => 'cnv_calling',
        pipeline    => $cnv_pipeline,
        datasource  => $fofn_ds,
        output_root => $output_dir,
        options     => {}
    );
    
    # figure out output files
    my (@merged_vcfs, @bychr_plot_files, @summary_files, @plot_files, @combine_files);
    foreach my $element_id (29 .. 31) {
        my @output_subdirs = output_subdirs($element_id, $cnv_setup->id);
        push(@merged_vcfs,      file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz'));
        push(@bychr_plot_files, file(@output_subdirs, '3_plot_polysomy',                             'copy_numbers.png'));
        push(@combine_files,    file(@output_subdirs, '5_combine_bcftools_cnvs',                     'combined_cnvs.txt'));
        if ($element_id == 30) {
            foreach my $sub_dir (qw(HPSI1013i-funy_2_QC1Hip-2197 HPSI1013i-funy_1_QC1Hip-2196 HPSI1013i-funy_3_QC1Hip-2198)) {
                my $this_dir = dir(@output_subdirs, '4_bcftools_cnv', $sub_dir);
                push(@summary_files, file($this_dir, 'summary.tab'));
                
                foreach my $chr (1 .. 22, 'X') {
                    next if "$chr" eq '6';
                    push(@plot_files, file($this_dir, "plot.HPSI1013pf-funy_QC1Hip-2195.$sub_dir.chr$chr.png"));
                }
            }
        }
    }
    ok handle_pipeline(@merged_vcfs, @bychr_plot_files, @summary_files, @plot_files, @combine_files), 'bcftools_cnv_caller pipeline ran ok and produced the expected output files';
    
    my $cnv_vrfile = VRPipe::File->get(path => $summary_files[0]);
    is $cnv_vrfile->meta_value('sample_control'), 'HPSI1013pf-funy_QC1Hip-2195', 'sample_control metadata exists on the file';
    
    my $vrtrack       = VRPipe::Schema->create("VRTrack");
    my $combined_node = $vrtrack->get_file($combine_files[0]);
    my @cnv_nodes     = $combined_node->related(outgoing => { namespace => 'VRTrack', label => 'CNVs' });
    is scalar(@cnv_nodes), 4, 'the combined cnvs file was parsed in to 4 CNVs nodes';
    my ($cnv_node) = grep { $_->md5_sample =~ /QC1Hip-2192/ } @cnv_nodes;
    my $data = $vrtrack->graph->json_decode($cnv_node->properties->{data});
    is_deeply $data, { ND => "0", LD => "0", SD => "0" }, 'one of the CNVs nodes had the expected data';
    my ($sample_node) = $cnv_node->related(incoming => { type => 'cnv_calls' });
    is $sample_node->name, 'QC1Hip-2192', 'the CNVs node was related to the correct sample';
    
    ok my $plot_node = $vrtrack->get_file($plot_files[0]), 'one of the cnv plots was in the graph db';
    ($sample_node) = $plot_node->related(incoming => { type => 'cnv_plot', namespace => 'VRTrack', label => 'Sample' });
    is $sample_node->name, 'QC1Hip-2197', 'the plot node was related to the correct sample';
    
    ($plot_node) = $sample_node->related(outgoing => { type => 'copy_number_by_chromosome_plot' });
    ok $plot_node, 'a copy_number_by_chromosome_plot node was associated with a sample';
}

# we'll also take the opportunity to test the loh caller pipeline, since that
# also uses files from the genome studio import
SKIP: {
    my $num_tests = 9;
    skip "hipsci loh calling tests disabled without hipsci_loh_caller.pl in your path", $num_tests unless can_execute('hipsci_loh_caller.pl');
    
    $output_dir = get_output_dir('hipsci_loh_caller');
    
    # check pipeline has correct steps
    ok my $loh_pipeline = VRPipe::Pipeline->create(name => 'hipsci_loh_caller'), 'able to get the hipsci_loh_caller pipeline';
    my @step_names;
    foreach my $stepmember ($loh_pipeline->step_members) {
        push(@step_names, $stepmember->step->name);
    }
    is_deeply \@step_names, [qw(vcf_merge_different_samples_control_aware hipsci_loh_caller)], 'the hipsci_loh_caller pipeline has the correct steps';
    
    my $loh_setup = VRPipe::PipelineSetup->create(
        name        => 'loh_calling',
        pipeline    => $loh_pipeline,
        datasource  => $fofn_ds,
        output_root => $output_dir,
        options     => {}
    );
    
    # figure out output files
    my (@merged_vcfs, @loh_files_with_results, @loh_files_no_results);
    foreach my $element_id (29 .. 31) {
        my @output_subdirs = output_subdirs($element_id, $loh_setup->id);
        push(@merged_vcfs, file(@output_subdirs, '1_vcf_merge_different_samples_control_aware', 'merged.vcf.gz'));
        
        my $result_file = file(@output_subdirs, '2_hipsci_loh_caller', 'merged.txt');
        if ($element_id == 31) {
            push(@loh_files_no_results, $result_file);
        }
        else {
            push(@loh_files_with_results, $result_file);
        }
    }
    ok handle_pipeline(@merged_vcfs, @loh_files_with_results), 'hipsci_loh_caller pipeline ran ok and produced the expected output files';
    
    my $created_empty = 0;
    foreach my $file (@loh_files_no_results) {
        $created_empty++ if (-e $file && !-s $file);
    }
    is $created_empty, 1, 'it also produced empty result files for the good cohorts';
    
    my $loh_vrfile = VRPipe::File->get(path => $loh_files_with_results[0]);
    is $loh_vrfile->meta_value('sample_control'), 'HPSI1013pf-garx_QC1Hip-2191', 'sample_control metadata exists on the file';
    
    my $vrtrack       = VRPipe::Schema->create("VRTrack");
    my $loh_file_node = $vrtrack->get_file($loh_files_with_results[0]);
    my @loh_nodes     = $loh_file_node->related(outgoing => { namespace => 'VRTrack', label => 'LOH' });
    is scalar(@loh_nodes), 1, 'the loh merged calls file was parsed in to 1 LOH node';
    my ($loh_node) = grep { $_->md5_sample =~ /QC1Hip-2192/ } @loh_nodes;
    my $data = $vrtrack->graph->json_decode($loh_node->properties->{data});
    is_deeply [$data->[0], $data->[-1], scalar(@$data)], [{ chr => 1, start => 10353112, end => 19765518, count => 15, type => 'SNPs' }, { chr => 22, start => 21928641, end => 37462926, count => 14, type => 'SNPs' }, 50], 'one of the LOH nodes had the expected data';
    my ($sample_node) = $loh_node->related(incoming => { type => 'loh_calls' });
    is $sample_node->name, 'QC1Hip-2192', 'the LOH node was related to the correct sample';
}

finish;
exit;

sub vcf_lines {
    my $vcf_file = shift;
    my @lines;
    open(my $fh, "zcat $vcf_file |");
    while (<$fh>) {
        next if /^##fileDate/;
        push(@lines, $_);
    }
    close($fh);
    return @lines;
}
