#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES GATK)],
        required_exe => [qw(samtools)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('bqsr_and_reduce');

ok my $bqsr_pipeline = VRPipe::Pipeline->create(name => 'bam_base_quality_score_recalibration_gatk_v2'), 'able to get the bam_base_quality_score_recalibration_gatk_v2 pipeline';
ok my $reduce_pipeline = VRPipe::Pipeline->create(name => 'bam_reduce_reads'), 'able to get the bam_reduce_reads pipeline';

my @s_names;
foreach my $stepmember ($bqsr_pipeline->steps, $reduce_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(
  bam_index
  gatk_base_recalibrator
  gatk_print_reads_with_bqsr
  bam_index
  bam_index
  gatk_reduce_reads
  bam_index);
is_deeply \@s_names, \@expected_step_names, 'the pipelines have the correct steps';

my $ref_fa_source   = file(qw(t data S_suis_P17.fa));
my $ref_dict_source = file(qw(t data S_suis_P17.fa.dict));
my $ref_dir         = dir($output_dir, 'ref');
$bqsr_pipeline->make_path($ref_dir);
my $ref_fa   = file($ref_dir, 'S_suis_P17.fa')->stringify;
my $ref_dict = file($ref_dir, 'S_suis_P17.dict')->stringify;
copy($ref_fa_source,   $ref_fa);
copy($ref_dict_source, $ref_dict);

my $res_dir = dir($output_dir, 'resources');
$bqsr_pipeline->make_path($res_dir);

my $known_sites_source = file(qw(t data known_sites.vcf.gz));
my $known_sites = file($res_dir, 'known_sites.vcf.gz')->stringify;
copy($known_sites_source,          $known_sites);
copy($known_sites_source . '.tbi', $known_sites . '.tbi');

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data improvement_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($output_dir, 'improvement_datasource.fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
while (<$ifh>) {
    chomp;
    my $source = file($_);
    my $dest = file($output_dir, $source->basename);
    copy($source, $dest);
    print $ofh $dest, "\n";
}
$orig_fofn_file->close;
$fofn_file->close;

my $bqsr_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 's_suis bqsr',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        source  => $fofn_file->path->stringify,
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $bqsr_pipeline,
    options     => {
        reference_fasta           => $ref_fa,
        base_recalibrator_options => "-l INFO -knownSites $known_sites -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate",
        cleanup                   => 0,
    }
);

my $reduce_pipelinesetup = VRPipe::PipelineSetup->create(
    name       => 's_suis reduce reads',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => '1[3]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $reduce_pipeline,
    options     => {
        reference_fasta => $ref_fa,
        cleanup         => 0,
    }
);

ok handle_pipeline(), 'pipelines ran okay';

my @output_files;
my @files = ('2822_7.pe.bam', '2822_6.pe.bam', '2822_6.se.bam', '2823_4.pe.bam', '8324_8.pe.bam');
my $element_id = 0;
foreach my $file (@files) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 1);
    $file =~ s/bam$/recal_data.grp/;
    push @output_files, file(@output_subdirs, '2_gatk_base_recalibrator', $file);
    $file =~ s/recal\_data\.grp$/recal.bam/;
    push @output_files, file(@output_subdirs, '3_gatk_print_reads_with_bqsr', $file);
    push @output_files, file(@output_subdirs, '3_gatk_print_reads_with_bqsr', $file . '.bai');
    
    foreach my $element (@{ get_elements($reduce_pipelinesetup->datasource) }) {
        my $path = shift @{ result_with_inflated_paths($element)->{paths} };
        next unless ($path eq $output_files[-2]);
        @output_subdirs = output_subdirs($element->id, 2);
        $file =~ s/bam$/reduced.bam/;
        push @output_files, file(@output_subdirs, '2_gatk_reduce_reads', $file);
        push @output_files, file(@output_subdirs, '2_gatk_reduce_reads', $file . '.bai');
    }
}
is scalar @output_files, 25, 'checking for the correct number of files';
ok handle_pipeline(@output_files), 'all expected output files were created';

finish;
