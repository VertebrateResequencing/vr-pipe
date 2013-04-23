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

my $output_dir = get_output_dir('convex_read_depth_generation_pipeline');

# convex read depths pipeline
#############################
ok my $pipeline1 = VRPipe::Pipeline->create(name => 'convex_read_depth_generation'), 'able to create the convex_read_depth_generation pipeline';

my @s_names;
foreach my $stepmember ($pipeline1->step_members) {
    push(@s_names, $stepmember->step->name);
}

my @expected_step_names = qw(convex_read_depth bam_metadata_with_sex);
is_deeply \@s_names, \@expected_step_names, 'the rd pipeline has the correct steps';

my $convex_r_libs = $ENV{CONVEX_R_LIB};
my $convex_home   = "$convex_r_libs/CoNVex";
my $classpath     = "$convex_home/java/lib/CoNVex.jar:$convex_home/java/lib/args4j-2.0.12.jar:$convex_home/java/lib/sam-1.67.jar";

my $regions_file = file(qw(t data cnv hs_chr20.convex.regions))->absolute->stringify;

my $fofn = file(qw(t data cnv cnv.bam.fofn))->absolute;

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

my (@output_files, @final_files);
my $element_id = 0;
my $f_fofn     = VRPipe::File->create(path => $fofn->absolute);
my $oh         = $f_fofn->openr;
while (<$oh>) {
    chomp;
    my $f  = file($_);
    my $fn = $f->basename;
    $fn =~ s/\.bam$/.rd.txt/;
    
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '1_convex_read_depth', $fn));
}
close($oh);

ok handle_pipeline(@output_files, @final_files), 'rd pipeline ran and created all expected output files';

# convex l2r pipeline
#####################
ok my $pipeline2 = VRPipe::Pipeline->create(name => 'convex_l2r_bp_generation'), 'able to create the convex_l2r_bp_generation pipeline';

@s_names = ();
foreach my $stepmember ($pipeline2->step_members) {
    push(@s_names, $stepmember->step->name);
}

@expected_step_names = qw(convex_breakpoints convex_L2R);
is_deeply \@s_names, \@expected_step_names, 'the l2r pipeline has the correct steps';

my $bp_file_name = "$output_dir/breakpoints.txt";

my $pipelinesetup2 = VRPipe::PipelineSetup->create(
    name       => 'convex_l2r_bp_generation_pipeline',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'convex_read_depth_generation_pipeline[convex_read_depth]',
        options => { metadata_keys => 'batch' }
    ),
    
    output_root => $output_dir,
    pipeline    => $pipeline2,
    options     => {
        cleanup             => 0,
        regions_file        => $regions_file,
        convex_rscript_path => "$convex_home/Rbatch",
        rscript_cmd         => '/software/bin/Rscript --vanilla',
        r_libs              => $convex_r_libs,
        includeChrX         => 0,
        bp_file_name        => $bp_file_name,
    }
);

for (my $i = 0; $i < @output_files; $i++) {
    $output_files[$i] =~ s/rd\.txt/l2r.txt/;
}
push @output_files, $bp_file_name;

ok handle_pipeline(@output_files, @final_files), 'l2r pipeline ran and created all expected output files';

# convex cnv call pipeline
#############################
ok my $pipeline3 = VRPipe::Pipeline->create(name => 'convex_cnv_calling'), 'able to create the convex_cnv_calling pipeline';

@s_names = ();
foreach my $stepmember ($pipeline3->step_members) {
    push(@s_names, $stepmember->step->name);
}

@expected_step_names = qw(convex_gam_correction convex_cnv_call);
is_deeply \@s_names, \@expected_step_names, 'the rd pipeline has the correct steps';

my $features_file       = file(qw(t data cnv hs_chr20.convex.features))->absolute->stringify;
my $centromere_reg_file = file(qw(t data cnv hs_chr20.convex.aps_table_hg19.txt))->absolute->stringify;

my $pipelinesetup3 = VRPipe::PipelineSetup->create(
    name        => 'convex_cnv_calling_pipeline',
    datasource  => VRPipe::DataSource->create(type => 'vrpipe', method => 'all', source => 'convex_read_depth_generation_pipeline[convex_read_depth]'),
    output_root => $output_dir,
    pipeline    => $pipeline3,
    options     => {
        cleanup             => 0,
        convex_rscript_path => "$convex_home/Rbatch",
        rscript_cmd         => '/software/bin/Rscript --vanilla',
        r_libs              => $convex_r_libs,
        features_file       => $features_file,
        breakpoints_file    => $bp_file_name,
        sw_exec             => "$convex_home/exec/swa_lin64",
        centromere_reg_file => $centromere_reg_file,
    }
);

my $i = 0;
my @cnv_output_files;
foreach my $de (@{ get_elements($pipelinesetup3->datasource) }) {
    my @output_subdirs = output_subdirs($de->id, $pipelinesetup3->id);
    
    my $basename = file($output_files[$i])->basename;
    $basename =~ s/l2r\.txt/gam.txt/;
    push(@cnv_output_files, file(@output_subdirs, '1_convex_gam_correction', $basename));
    $basename =~ s/gam\.txt/cnv_calls.txt/;
    push(@cnv_output_files, file(@output_subdirs, '2_convex_cnv_call', $basename));
    
    $i++;
}

ok handle_pipeline(@cnv_output_files, @final_files), 'cnv call pipeline ran and created all expected output files';

done_testing;
