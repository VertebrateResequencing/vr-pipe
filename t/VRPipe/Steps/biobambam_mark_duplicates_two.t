#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 2;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::biobambam_mark_duplicates_two');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('biobambam_mark_duplicates_two', 'cram_files');
is_deeply [$step->id, $step->description], [1, 'Marks duplicates in CRAM files using biobambam'], 'biobambam_mark_duplicates_two step created and has correct description';

my $setup = VRPipe::PipelineSetup->create(
    name       => 'biobambam_mark_duplicates_two',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'grouped_by_metadata',
        source  => file(qw(t data cram.fofn_with_metadata))->absolute->stringify,
        options => {
            metadata_keys = 'sample',
        }
    },
    output_root => $output_dir,
    pipeline    => $pipeline,
    options     => {}
};

ok handle_pipeline(), 'biobambam_mark_duplicates_two pipeline ran ok';

finish;
