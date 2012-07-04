#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 3;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES)],
		    required_exe => [qw(smalt)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::smalt_index');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('smalt_index', '');
is_deeply [$step->id, $step->description], [1, 'Indexes a reference genome fasta file, making it suitable for use in subsequent smalt mapping'], 'smalt_index step created and has correct description';

# test using the class methods directly
my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($output_dir, 'ref');
$pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->create(name => 'smalt_setup',
                                       datasource => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data decompress_datasource.fofn))->absolute),
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => { reference_fasta => $ref_fa, smalt_index_options => '-k 13 -s 4'});

my @outfiles = (file($ref_dir, 'S_suis_P17.fa_k13s4.smi'), file($ref_dir, 'S_suis_P17.fa_k13s4.sma'));
ok handle_pipeline(@outfiles), 'single-step pipeline ran ok';

finish;
