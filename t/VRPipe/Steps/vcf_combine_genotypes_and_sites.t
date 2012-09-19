#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(tabix)],
        required_module => [qw(Vcf)]
    );
    use TestPipelines;
    
    use_ok('VRPipe::Steps::vcf_combine_genotypes_and_sites');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('vcf_combine_genotypes_and_sites', 'genotype_vcf_files');
is_deeply [$step->id, $step->description], [1, 'Merge a sites only VCF file with the genotypes found in another VCF file, keeping the FILTER, ID and INFO fields found in the sites only file'], 'vcf_combine_genotypes_and_sites step created and has correct description';

VRPipe::PipelineSetup->create(
    name        => 'vcf_combine_genotypes_and_sites step test',
    datasource  => VRPipe::DataSource->create(type => 'fofn_with_metadata', method => 'all', source => file(qw(t data combine_site_and_genotypes_datasource.fofn))->absolute),
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {
        sites_vcf_path => file(qw(t data sites.vcf.gz))->absolute->stringify,
    }
);

my @output_files;
my @output_subdirs = output_subdirs(1, 1);
foreach my $eid (1, 2) {
    my @output_subdirs = output_subdirs($eid, 1);
    push(@output_files, file(@output_subdirs, '1_vcf_combine_genotypes_and_sites', 'combined.vcf.gz'));
}
ok handle_pipeline(@output_files), 'single step vcf_combine_genotypes_and_sites pipeline ran and created all expected output files';

is num_samples($output_files[0]->absolute->stringify), 16, 'correct number of samples in first output file';

sub num_samples {
    my $vcf_path = shift;
    my $vcf = Vcf->new(file => $vcf_path);
    $vcf->parse_header();
    my $samples = scalar @{ $vcf->{columns} } - 9;
    $vcf->close;
    return $samples;
}

finish;
