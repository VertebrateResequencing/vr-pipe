#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(bcftools tabix bgzip)]
    );
    use TestPipelines;
}

my $ref_dir      = get_output_dir('reference');
my $ref_vcf_orig = file(qw(t data chr20.1kg_phase1.GBR.vcf.gz));
my $ref_vcf      = file($ref_dir, 'chr20.1kg_phase1.GBR.vcf.gz')->stringify;
copy($ref_vcf_orig,       $ref_vcf);
copy("$ref_vcf_orig.tbi", "$ref_vcf.tbi");

SKIP: {
    my $num_tests = 3;
    skip "shapeit test disabled without shapeit in your path", $num_tests unless can_execute('shapeit');
    
    my $output_dir = get_output_dir('shapeit');
    # check pipeline has correct steps
    ok my $shapeit_pipeline = VRPipe::Pipeline->create(name => 'shapeit'), 'able to get the shapeit pipeline';
    my @sb_names;
    foreach my $stepmember ($shapeit_pipeline->step_members) {
        push(@sb_names, $stepmember->step->name);
    }
    
    is_deeply \@sb_names, [qw(define_vcf_chunks create_imputation_ref_panel shapeit bcftools_concat)], 'the shapeit pipeline has the correct steps';
    
    my $shapeit_setup = VRPipe::PipelineSetup->create(
        name       => 'shapeit setup',
        pipeline   => $shapeit_pipeline,
        datasource => VRPipe::DataSource->create(
            type   => 'fofn',
            method => 'all',
            source => file(qw(t data imputation.vcf.fofn)),
        ),
        options => {
            chunk_nsites         => 1000,
            buffer_nsites        => 200,
            shapeit_options      => '--seed 123456789',
            merge_by_chromosome  => 1,
            regions              => '20:797450-1000000,1:61284739-249250621,20:128473-192506',
            bcftools_concat_opts => '--ligate --allow-overlaps'
        },
        output_root => $output_dir
    );
    
    my @final_files;
    foreach my $element (@{ get_elements($shapeit_setup->datasource) }) {
        my @output_dirs = output_subdirs($element->id, $shapeit_setup->id);
        push(@final_files, file(@output_dirs, '4_bcftools_concat', 'merged.chr20.vcf.gz.tbi'));
    }
    ok handle_pipeline(@final_files), 'shapeit pipeline ran ok';

}

SKIP: {
    my $num_tests = 3;
    skip "impute2 test disabled without impute2 in your path", $num_tests unless can_execute('impute2');
    
    my $output_dir = get_output_dir('impute2');
    # check pipeline has correct steps
    ok my $impute2_pipeline = VRPipe::Pipeline->create(name => 'impute2'), 'able to get the impute2 pipeline';
    my @sb_names;
    foreach my $stepmember ($impute2_pipeline->step_members) {
        push(@sb_names, $stepmember->step->name);
    }
    
    is_deeply \@sb_names, [qw(define_vcf_chunks create_imputation_ref_panel impute2 bcftools_concat)], 'the impute2 pipeline has the correct steps';
    
    my $impute2_setup = VRPipe::PipelineSetup->create(
        name       => 'impute2 setup',
        pipeline   => $impute2_pipeline,
        datasource => VRPipe::DataSource->create(
            type   => 'fofn',
            method => 'all',
            source => file(qw(t data imputation.vcf.fofn)),
        ),
        options => {
            ref_vcf       => $ref_vcf,
            chunk_by_ref  => 1,
            chunk_nsites  => 1000,
            buffer_nsites => 200
        },
        output_root => $output_dir
    );
    
    my @final_files;
    foreach my $element (@{ get_elements($impute2_setup->datasource) }) {
        my @output_dirs = output_subdirs($element->id, $impute2_setup->id);
        push(@final_files, file(@output_dirs, '4_bcftools_concat', 'merged.vcf.gz.tbi'));
    }
    ok handle_pipeline(@final_files), 'impute2 pipeline ran ok';

}

finish;
exit;
