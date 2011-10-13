#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 5;
    
    use_ok('VRPipe::Persistent::Schema');
    use_ok('VRPipe::Steps::bam_strip_tags');
    
    use TestPipelines;
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_strip_tags', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Strips tags from bam files'], 'bam_strip_tags step created and has correct description';

# test using the class methods directly
my $test_bam = file(qw(t data 2822_6.pe.bam))->absolute;

my @tags = qw(OQ XM XG XO E2);
my $tag_bam = file(qw(t data 2822_6.pe.bam))->absolute;
my $strip_bam = file($output_dir, 'strip.bam');
my $ok = VRPipe::Steps::bam_strip_tags->tag_strip($tag_bam, $strip_bam, tags_to_strip => [@tags]);
is $ok, 1,  'tag_strip() ran okay';

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->get(name => 'strip_setup',
                                       datasource => VRPipe::DataSource->get(type => 'fofn', method => 'all', source => file(qw(t data improvement_datasource.fofn))->absolute),
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => { bam_tags_to_strip => [qw(OQ XM XG XO E2)] });

ok handle_pipeline(), 'single-step pipeline ran ok';

finish;
