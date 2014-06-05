#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 9;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES CONVEX_R_LIB)],
        max_retries  => 1,
        required_exe => [qw(samtools R)]
    );
    use TestPipelines;
}

my $output_dir = get_output_dir('convex_pipeline');

# convex read depths pipeline
#############################
ok my $pipeline1 = VRPipe::Pipeline->create(name => 'convex_read_depth_generation'), 'able to create the convex_read_depth_generation pipeline';

my @s_names;
foreach my $stepmember ($pipeline1->step_members) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(bam_metadata_with_sex convex_read_depth);
is_deeply \@s_names, \@expected_step_names, 'the convex_read_depth_generation pipeline has the correct steps';

my $convex_r_libs = $ENV{CONVEX_R_LIB};
my $convex_home   = "$convex_r_libs/CoNVex";
my $classpath     = "$convex_home/java/lib/CoNVex.jar:$convex_home/java/lib/args4j-2.0.12.jar:$convex_home/java/lib/sam-1.67.jar";

my $regions_file = file(qw(t data cnv hs_chr20.convex.regions))->absolute->stringify;

my $fofn = file(qw(t data cnv cnv.bam.fofn));

my $pipelinesetup1 = VRPipe::PipelineSetup->create(
    name        => 'convex_read_depth_generation_pipeline',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => $fofn),
    output_root => $output_dir,
    pipeline    => $pipeline1,
    options     => {
        cleanup          => 0,
        convex_classpath => $classpath,
        regions_file     => $regions_file,
        assumed_sex      => 'U',
    }
);

my @output_files;
my $element_id = 0;
my $f_fofn     = VRPipe::File->create(path => $fofn->absolute);
my $oh         = $f_fofn->openr;
while (<$oh>) {
    chomp;
    my $f  = file($_);
    my $fn = $f->basename;
    $fn =~ s/\.bam$/.rd/;
    
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '2_convex_read_depth', $fn));
}
close($oh);

ok handle_pipeline(@output_files), 'convex_read_depth_generation ran and created all expected output files';

# convex l2r pipeline
#####################
ok my $pipeline2 = VRPipe::Pipeline->create(name => 'convex_l2r_bp_generation'), 'able to create the convex_l2r_bp_generation pipeline';

@s_names = ();
foreach my $stepmember ($pipeline2->step_members) {
    push(@s_names, $stepmember->step->name);
}

@expected_step_names = qw(convex_breakpoints convex_L2R);
is_deeply \@s_names, \@expected_step_names, 'the convex_l2r_bp_generation pipeline has the correct steps';

my $pipelinesetup2 = VRPipe::PipelineSetup->create(
    name       => 'convex_l2r_bp_generation_pipeline',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_all',
        source  => 'convex_read_depth_generation_pipeline[0,2]',
        options => {}
    ),
    
    output_root => $output_dir,
    pipeline    => $pipeline2,
    options     => {
        cleanup             => 0,
        regions_file        => $regions_file,
        convex_rscript_path => "$convex_home/Rbatch",
        rscript_cmd         => 'Rscript --vanilla',
        r_libs              => $convex_r_libs,
        includeChrX         => 0,
        minSamples          => 2
    }
);

for (my $i = 0; $i < @output_files; $i++) {
    $output_files[$i] =~ s/\.rd/.l2r/;
}

my $de = ${ get_elements($pipelinesetup2->datasource) }[0];
my @output_subdirs = output_subdirs($de->id, $pipelinesetup2->id);
push(@output_files, file(@output_subdirs, '2_convex_L2R', 'sample_info.txt')->absolute->stringify);
my $features_file     = file(@output_subdirs, '2_convex_L2R', 'features.fts')->absolute->stringify;
my $corr_matrix_file  = file(@output_subdirs, '2_convex_L2R', 'corr_matrix.corr')->absolute->stringify;
my $sample_means_file = file(@output_subdirs, '2_convex_L2R', 'sample_means.savg')->absolute->stringify;
push(@output_files, $features_file, $corr_matrix_file, $sample_means_file);

ok handle_pipeline(@output_files), 'convex_l2r_bp_generation pipeline ran and created all expected output files';

# convex cnv call pipeline
#############################
ok my $pipeline3 = VRPipe::Pipeline->create(name => 'convex_cnv_calling'), 'able to create the convex_cnv_calling pipeline';

@s_names = ();
foreach my $stepmember ($pipeline3->step_members) {
    push(@s_names, $stepmember->step->name);
}

@expected_step_names = qw(convex_gam_correction convex_cnv_call);
is_deeply \@s_names, \@expected_step_names, 'the convex_cnv_calling pipeline has the correct steps';

my $centromere_reg_file = file(qw(t data cnv hs_chr20.convex.aps_table_hg19.txt))->absolute->stringify;

my $pipelinesetup3 = VRPipe::PipelineSetup->create(
    name       => 'convex_cnv_calling_pipeline',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'convex_read_depth_generation_pipeline[convex_read_depth]|convex_l2r_bp_generation_pipeline[convex_L2R:l2r_files]',
        options => {
            metadata_keys           => 'sample',
            include_in_all_elements => 'convex_l2r_bp_generation_pipeline[convex_breakpoints:breakpoints_file]|convex_l2r_bp_generation_pipeline[convex_L2R:features_file]'
        }
    ),
    output_root => $output_dir,
    pipeline    => $pipeline3,
    options     => {
        cleanup             => 0,
        convex_rscript_path => "$convex_home/Rbatch",
        rscript_cmd         => 'Rscript --vanilla',
        r_libs              => $convex_r_libs,
        sw_exec             => "$convex_home/exec/swa_lin64",
        centromere_reg_file => $centromere_reg_file,
    }
);

my $i = 0;
my @cnv_output_files;
foreach my $de (@{ get_elements($pipelinesetup3->datasource) }) {
    @output_subdirs = output_subdirs($de->id, $pipelinesetup3->id);
    
    my $basename = file($output_files[$i])->basename;
    $basename =~ s/l2r$/gam/;
    push(@cnv_output_files, file(@output_subdirs, '1_convex_gam_correction', $basename));
    $basename =~ s/gam$/cnv/;
    push(@cnv_output_files, file(@output_subdirs, '2_convex_cnv_call', $basename));
    
    $i++;
}

ok handle_pipeline(@cnv_output_files), 'cnv call pipeline ran and created all expected output files';

# convex_cnv_postprocess pipeline
#############################
# ok my $pipeline4 = VRPipe::Pipeline->create(name => 'convex_cnv_postprocess'), 'able to create the convex_cnv_postprocess pipeline';

# @s_names = ();
# foreach my $stepmember ($pipeline4->step_members) {
#     push(@s_names, $stepmember->step->name);
# }

# @expected_step_names = qw(convex_calls_with_means_mads convex_plots);
# is_deeply \@s_names, \@expected_step_names, 'the convex_cnv_postprocess pipeline has the correct steps';

# my $pipelinesetup4 = VRPipe::PipelineSetup->create(
#     name        => 'convex_postprocess_pipeline',
#     datasource  => VRPipe::DataSource->create(
#         type => 'vrpipe',
#         method => 'group_all',
#         source => '1[0,2]|2[2:l2r_files,2:corr_matrix_file,2:features_file]|3[1,2]', # bam,rd,l2r,gam,cnv
#     ),
#     output_root => $output_dir,
#     pipeline    => $pipeline4,
#     options     => {
#         cleanup             => 0,
#         convex_rscript_path => "$convex_home/Rbatch",
#         convex_classpath    => $classpath,
#         rscript_cmd         => 'Rscript --vanilla',
#         r_libs              => $convex_r_libs,
#         regions_file        => $regions_file,
#     }
# );

# $i = 0;
# my @mm_output_files;
# foreach my $cnv_file (@cnv_output_files) {
#     my $mm_file = $cnv_file;
#     $mm_file =~ s/cnv$/mm.cnv/;
#     push @mm_output_files, $mm_file;
# }

# ok handle_pipeline(@mm_output_files), 'convex_cnv_postprocess pipeline ran and created all expected output files';

# my @plot_files;
# foreach my $de (@{ get_elements($pipelinesetup4->datasource) }) {
#     @output_subdirs = output_subdirs($de->id, $pipelinesetup4->id);
#     push(@plot_files, file(@output_subdirs, '2_convex_plots', 'CNVstats_CallsperSample.png'));
#     push(@plot_files, file(@output_subdirs, '2_convex_plots', 'CNVstats_DelDupRatio.png'));
# }

# ok handle_pipeline(@plot_files), 'plot files created as expected';

done_testing;
