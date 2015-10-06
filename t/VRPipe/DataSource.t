#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use Parallel::ForkManager;

BEGIN {
    use Test::Most tests => 183;
    use VRPipeTest;
    use TestPipelines;
    
    use_ok('VRPipe::DataSourceFactory');
    use_ok('VRPipe::Schema');
}

# list
ok my $ds = VRPipe::DataSource->create(
    type    => 'list',
    method  => 'all',
    source  => file(qw(t data datasource.list))->absolute->stringify,
    options => {}
  ),
  'could create a list datasource';

my @results;
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ line => 'foo' }, { line => 'bar' }, { line => 'henry' }], 'got all results correctly';

$ds = VRPipe::DataSource->create(
    type    => 'list',
    method  => 'all',
    source  => file(qw(t data datasource.list))->absolute->stringify,
    options => { skip_comments => 0 }
);

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ line => 'foo' }, { line => 'bar' }, { line => '# comment' }, { line => 'henry' }], 'got even more results with extra options';

# fofn
ok $ds = VRPipe::DataSource->create(
    type    => 'fofn',
    method  => 'all',
    source  => file(qw(t data datasource.fofn))->absolute->stringify,
    options => {}
  ),
  'could create a fofn datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file('t', 'data', 'file.bam')->absolute] }, { paths => [file('t', 'data', 'file.cat')->absolute] }, { paths => [file('t', 'data', 'file.txt')->absolute] }], 'got correct results for fofn all';

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn',
    method  => 'group_all',
    source  => file(qw(t data datasource.fofn))->absolute->stringify,
    options => {}
  ),
  'could create a fofn datasource with group_all method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file('t', 'data', 'file.bam')->absolute, file('t', 'data', 'file.cat')->absolute, file('t', 'data', 'file.txt')->absolute] }], 'got correct results for fofn group_all';

# test that _prepare_elements_and_states only happens once simultaneously
my $fm           = Parallel::ForkManager->new(4);
my $num_prepared = 0;
$fm->run_on_finish(
    sub {
        my ($pid, $prepared) = @_;
        $num_prepared++ if $prepared;
    }
);
for (1 .. 4) {
    $fm->start and next;
    my $child_ds = VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        source  => file(qw(t data datasource.fofn2))->absolute->stringify,
        options => {}
    );
    # force them to be considered changed, as if this was a ds that took a long
    # time to prepare the first time
    $child_ds->_changed_marker('foo');
    $child_ds->_source_instance->_changed_marker('foo');
    my $prepared = $child_ds->_prepare_elements_and_states;
    $fm->finish($prepared);
}
$fm->wait_all_children;
is $num_prepared, 1, '_prepare_elements_and_states only occurs once at a time';

# properly test that datasource updates work correctly when we use the same
# datasource instance
{
    # create a simple fofn ds with multiple elements to test with
    my $output_root       = get_output_dir('datasource_updating_fofn_output');
    my $source_file       = file($output_root, 'source.fofn');
    my $source_inputs_dir = dir($output_root, 'inputs');
    mkdir($source_inputs_dir);
    my $sfh          = $source_file->openw;
    my $num_elements = 5;
    for my $i (1 .. $num_elements) {
        my $input_file = file($source_inputs_dir, 'input.' . $i . '.txt');
        my $ifh = $input_file->openw;
        print $ifh "$i\n";
        close($ifh);
        print $sfh "$input_file\n";
    }
    close($sfh);
    my $ds = VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        source  => $source_file->stringify,
        options => {}
    );
    
    my $element_count = @{ get_elements($ds) };
    is $element_count, $num_elements, "our test datasource starts with $num_elements dataelements";
    
    # now update the ds to have more elements
    open($sfh, '>>', $source_file);
    my $new_num_elements = ($num_elements / 5) + $num_elements;
    for my $i (($num_elements + 1) .. $new_num_elements) {
        my $input_file = file($source_inputs_dir, 'input.' . $i . '.txt');
        my $ifh = $input_file->openw;
        print $ifh "$i\n";
        close($ifh);
        print $sfh "$input_file\n";
    }
    close($sfh);
    
    # check the element count is correct
    $element_count = @{ get_elements($ds) };
    is $element_count, $new_num_elements, "the test datasource correctly updated to $new_num_elements dataelements even when we reused the same instance";
}

# delimited
ok $ds = VRPipe::DataSource->create(
    type    => 'delimited',
    method  => 'grouped_single_column',
    source  => file(qw(t data datasource.fastqs))->absolute->stringify,
    options => {
        delimiter => "\t",
        group_by  => 1,
        column    => 2
    }
  ),
  'could create a delimited datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file(qw(t data 2822_6_1.fastq))->absolute, file(qw(t data 2822_6_2.fastq))->absolute], group => '2822_6' }, { paths => [file(qw(t data 2822_7_1.fastq))->absolute, file(qw(t data 2822_7_2.fastq))->absolute], group => '2822_7' }, { paths => [file(qw(t data 2823_4_1.fastq))->absolute, file(qw(t data 2823_4_2.fastq))->absolute], group => '2823_4' }, { paths => [file(qw(t data 8324_8_1.fastq))->absolute, file(qw(t data 8324_8_2.fastq))->absolute], group => '8324_8' }], 'got correct results for delimited grouped_single_column';

# delimited all columns
ok $ds = VRPipe::DataSource->create(
    type    => 'delimited',
    method  => 'all_columns',
    source  => file(qw(t data datasource.2col))->absolute->stringify,
    options => { delimiter => "\t" }
  ),
  'could create a delimited datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file(qw(t data file.txt))->absolute, file(qw(t data file2.txt))->absolute] }, { paths => [file(qw(t data file3.txt))->absolute, file(qw(t data file.cat))->absolute] }], 'got correct results for delimited all_columns';

# test get_methods and method_options
$ds = $ds->_source_instance;
is_deeply [sort $ds->get_methods], [qw(all all_columns grouped_single_column single_column)], 'get_methods returned all the expected methods';
is_deeply [$ds->method_options('single_column')], [['named', 'delimiter', 1, undef, 'Str'], ['named', 'column', 1, undef, 'PositiveInt'], ['named', 'column_is_path', 0, 1, 'Bool']], 'method_options showed us what options the single_column method takes';

# fofn_with_metadata
my @fwm_paths = ('/a/path/7816_3#95.bam', '/a/path/7413_5#95.bam', '/a/path/8312_5#95.bam');
my %fwm_common_meta = (center_name => 'SC', study => 'ERP000979', platform => 'ILLUMINA');
VRPipe::File->create(path => $fwm_paths[0])->add_metadata({ library => 'foo' });
ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'all',
    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
    options => {}
  ),
  'could create a fofn_with_metadata datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply [@results, VRPipe::File->get(path => $fwm_paths[0])->metadata, VRPipe::File->get(path => $fwm_paths[1])->metadata, VRPipe::File->get(path => $fwm_paths[2])->metadata], [{ paths => [$fwm_paths[0]] }, { paths => [$fwm_paths[1]] }, { paths => [$fwm_paths[2]] }, { %fwm_common_meta, sample => 'JB953', library => '4858080', lane => '7816_3#95' }, { %fwm_common_meta, sample => 'JB953', library => '4074406', lane => '7413_5#95' }, { %fwm_common_meta, sample => 'JB951', library => '4074399', lane => '8312_5#95' }], 'got correct results for fofn_with_metadata all, and the metadata on the files was correct';

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'all',
    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
    options => { filter => 'sample#JB953' }
  ),
  'could create a fofn_with_metadata datasource and filter';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply [@results, VRPipe::File->get(path => $fwm_paths[0])->metadata, VRPipe::File->get(path => $fwm_paths[1])->metadata], [{ paths => [$fwm_paths[0]] }, { paths => [$fwm_paths[1]] }, { %fwm_common_meta, sample => 'JB953', library => '4858080', lane => '7816_3#95' }, { %fwm_common_meta, sample => 'JB953', library => '4074406', lane => '7413_5#95' }], 'got correct results for fofn_with_metadata all + filtering, and the metadata on the files was correct';

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'group_all',
    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
    options => {}
  ),
  'could create a fofn_with_metadata datasource with group_all method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply [@results, VRPipe::File->get(path => $fwm_paths[0])->metadata, VRPipe::File->get(path => $fwm_paths[1])->metadata, VRPipe::File->get(path => $fwm_paths[2])->metadata], [{ paths => [$fwm_paths[0], $fwm_paths[1], $fwm_paths[2]] }, { %fwm_common_meta, sample => 'JB953', library => '4858080', lane => '7816_3#95' }, { %fwm_common_meta, sample => 'JB953', library => '4074406', lane => '7413_5#95' }, { %fwm_common_meta, sample => 'JB951', library => '4074399', lane => '8312_5#95' }], 'got correct results for fofn_with_metadata group_all, and the metadata on the files was correct';

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'group_all',
    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
    options => { filter => 'sample#JB953' }
  ),
  'could create a fofn_with_metadata datasource with group_all method and filter';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply [@results, VRPipe::File->get(path => $fwm_paths[0])->metadata, VRPipe::File->get(path => $fwm_paths[1])->metadata], [{ paths => [$fwm_paths[0], $fwm_paths[1]] }, { %fwm_common_meta, sample => 'JB953', library => '4858080', lane => '7816_3#95' }, { %fwm_common_meta, sample => 'JB953', library => '4074406', lane => '7413_5#95' }], 'got correct results for fofn_with_metadata group_all + filtering, and the metadata on the files was correct';

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'grouped_by_metadata',
    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
    options => { metadata_keys => 'study|sample' }
  ),
  'could create a fofn_with_metadata grouped_by_metadata datasource';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply [sort { $a->{group} cmp $b->{group} } @results], [{ paths => [$fwm_paths[2]], group => 'ERP000979|JB951' }, { paths => [$fwm_paths[0], $fwm_paths[1]], group => 'ERP000979|JB953' }], 'got correct results for fofn_with_metadata grouped_by_metadata';

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_metadata',
    method  => 'grouped_by_metadata',
    source  => file(qw(t data datasource.fofn_with_metadata))->absolute->stringify,
    options => { metadata_keys => 'study|sample', filter => 'sample#JB953' }
  ),
  'could create a fofn_with_metadata grouped_by_metadata datasource and filter';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply [sort { $a->{group} cmp $b->{group} } @results], [{ paths => [$fwm_paths[0], $fwm_paths[1]], group => 'ERP000979|JB953' }], 'got correct results for fofn_with_metadata grouped_by_metadata + filtering';

# sequence_index
ok $ds = VRPipe::DataSource->create(
    type    => 'sequence_index',
    method  => 'lane_fastqs',
    source  => file(qw(t data datasource.sequence_index))->absolute->stringify,
    options => { local_root_dir => dir('./')->absolute->stringify }
  ),
  'could create a sequence_index datasource with lane_fastqs method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file(qw(t data 2822_6.fastq))->absolute, file(qw(t data 2822_6_1.fastq))->absolute, file(qw(t data 2822_6_2.fastq))->absolute], lane => '2822_6' }, { paths => [file(qw(t data 2822_7_1.fastq))->absolute, file(qw(t data 2822_7_2.fastq))->absolute], lane => '2822_7' }, { paths => [file(qw(t data 2823_4_1.fastq))->absolute, file(qw(t data 2823_4_2.fastq))->absolute], lane => '2823_4' }, { paths => [file(qw(t data 8324_8_1.fastq))->absolute, file(qw(t data 8324_8_2.fastq))->absolute], lane => '8324_8' }], 'got correct results for sequence_index lane_fastqs';
my $vrfile = VRPipe::File->get(path => file(qw(t data 2822_6_1.fastq))->absolute);
my $meta = $vrfile->metadata;
is_deeply $meta,
  {
    expected_md5   => 'f1826489facca0d0bdf02d9586b493f6',
    lane           => '2822_6',
    study          => 'STUDY01',
    study_name     => 'my study name',
    center_name    => 'SC',
    sample_id      => 'SAMPLEID01',
    sample         => 'SAMPLE01',
    population     => 'POP',
    platform       => 'ILLUMINA',
    library        => 'LIB01',
    insert_size    => 200,
    withdrawn      => 0,
    reads          => 200,
    bases          => 12200,
    analysis_group => 'low coverage',
    paired         => 1,
    mate           => file(qw(t data 2822_6_2.fastq))->absolute->stringify
  },
  'a VRPipe::File created by source has the correct metadata';

ok $ds = VRPipe::DataSource->create(
    type    => 'sequence_index',
    method  => 'sample_fastqs',
    source  => file(qw(t data datasource.sequence_index))->absolute->stringify,
    options => { local_root_dir => dir('./')->absolute->stringify }
  ),
  'could create a sequence_index datasource with sample_fastqs method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
is_deeply \@results, [{ paths => [file(qw(t data 2822_6.fastq))->absolute, file(qw(t data 2822_6_1.fastq))->absolute, file(qw(t data 2822_6_2.fastq))->absolute, file(qw(t data 2822_7_1.fastq))->absolute, file(qw(t data 2822_7_2.fastq))->absolute, file(qw(t data 2823_4_1.fastq))->absolute, file(qw(t data 2823_4_2.fastq))->absolute], sample => 'SAMPLE01' }, { paths => [file(qw(t data 8324_8_1.fastq))->absolute, file(qw(t data 8324_8_2.fastq))->absolute], sample => 'SAMPLE02' }], 'got correct results for sequence_index sample_fastqs';

my $fai      = file(qw(t data human_g1k_v37.fasta.fai))->absolute->stringify;
my $override = file(qw(t data wgs_calling_override_options))->absolute->stringify;

# genome chunking
my $chunks = [{ chrom => 11, from => 1, to => 10000000, seq_no => 1, chunk_override_file => $override }, { chrom => 11, from => 10000001, to => 20000000, seq_no => 2, chunk_override_file => $override }, { chrom => 11, from => 20000001, to => 30000000, seq_no => 3, chunk_override_file => $override }, { chrom => 11, from => 30000001, to => 40000000, seq_no => 4, chunk_override_file => $override }, { chrom => 11, from => 40000001, to => 50000000, seq_no => 5, chunk_override_file => $override }, { chrom => 11, from => 50000001, to => 60000000, seq_no => 6, chunk_override_file => $override }, { chrom => 11, from => 60000001, to => 70000000, seq_no => 7, chunk_override_file => $override }, { chrom => 11, from => 70000001, to => 80000000, seq_no => 8, chunk_override_file => $override }, { chrom => 11, from => 80000001, to => 90000000, seq_no => 9, chunk_override_file => $override }, { chrom => 11, from => 90000001, to => 100000000, seq_no => 10, chunk_override_file => $override }, { chrom => 11, from => 100000001, to => 110000000, seq_no => 11, chunk_override_file => $override }, { chrom => 11, from => 110000001, to => 120000000, seq_no => 12, chunk_override_file => $override }, { chrom => 11, from => 120000001, to => 130000000, seq_no => 13, chunk_override_file => $override }, { chrom => 11, from => 130000001, to => 135006516, seq_no => 14, chunk_override_file => $override }, { chrom => 20, from => 1, to => 10000000, seq_no => 15, chunk_override_file => $override }, { chrom => 20, from => 10000001, to => 20000000, seq_no => 16, chunk_override_file => $override }, { chrom => 20, from => 20000001, to => 30000000, seq_no => 17, chunk_override_file => $override }, { chrom => 20, from => 30000001, to => 40000000, seq_no => 18, chunk_override_file => $override }, { chrom => 20, from => 40000001, to => 50000000, seq_no => 19, chunk_override_file => $override }, { chrom => 20, from => 50000001, to => 60000000, seq_no => 20, chunk_override_file => $override }, { chrom => 20, from => 60000001, to => 63025520, seq_no => 21, chunk_override_file => $override }];

my $chunking_ds_options = { reference_index => $fai, chunk_override_file => $override, chrom_list => '11 20', chunk_size => 10000000 };
ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_genome_chunking',
    method  => 'group_all',
    source  => file(qw(t data bams.fofn))->absolute->stringify,
    options => $chunking_ds_options
  ),
  'could create a fofn_with_genome_chunking datasource with group_all method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
my @expected = ();
foreach my $chunk (@$chunks) {
    push(@expected, { paths => [file('t', 'data', 'NA19334.bam')->absolute, file('t', 'data', 'NA19381.bam')->absolute, file('t', 'data', 'NA20281.bam')->absolute], %$chunk });
}
is_deeply \@results, \@expected, 'got correct results for fofn_with_genome_chunking group_all method';

my $ds_si = $ds->_source_instance;
is_deeply [$ds_si->method_options('group_all')], [['named', 'reference_index', 1, undef, 'Str|File'], ['named', 'chunk_override_file', 1, undef, 'Str|File'], ['named', 'chunk_size', 1, '1000000', 'Int'], ['named', 'chunk_overlap', 1, '0', 'Int'], ['named', 'chrom_list', 0, undef, 'Str'], ['named', 'ploidy', 0, undef, 'Str|File'], ['named', 'target_regions', 0, undef, 'Str|File']], 'method_options call for fofn_with_genome_chunking datasource got correct result';

is $ds_si->method_description('group_all'), q[All files in the file will be grouped into a single element. Each dataelement will be duplicated in chunks across the genome. The option 'reference_index' is the absolute path to the fasta index (.fai) file associated with the reference fasta file, 'chunk_override_file' is a file defining chunk specific options that may be overridden (required, but may point to an empty file), 'chunk_size' the size of the chunks in bp, 'chunk_overlap' defines how much overlap to have beteen chunks, 'chrom_list' (a space separated list) will restrict to specified the chromosomes (must match chromosome names in dict file), 'ploidy' is an optional file specifying the ploidy to be used for males and females in defined regions of the genome, eg {default=>2, X=>[{ from=>1, to=>60_000, M=>1 },{ from=>2_699_521, to=>154_931_043, M=>1 },],Y=>[{ from=>1, to=>59_373_566, M=>1, F=>0 }]}.], 'method description for fofn_with_genome_chunking group_all method is correct';

my $chroms = [{ chrom => 11, from => 1, to => 135006516, seq_no => 1, chunk_override_file => $override }, { chrom => 20, from => 1, to => 63025520, seq_no => 2, chunk_override_file => $override }];

ok my $ds2 = VRPipe::DataSource->create(
    type    => 'fofn_with_genome_chunking',
    method  => 'group_all',
    source  => file(qw(t data bams.fofn))->absolute->stringify,
    options => { %$chunking_ds_options, chunk_size => 0 }
  ),
  'could create a fofn_with_genome_chunking datasource with group_all method and chunk_size 0';

@results = ();
foreach my $element (@{ get_elements($ds2) }) {
    push(@results, result_with_inflated_paths($element));
}
@expected = ();
foreach my $chrom (@$chroms) {
    push(@expected, { paths => [file('t', 'data', 'NA19334.bam')->absolute, file('t', 'data', 'NA19381.bam')->absolute, file('t', 'data', 'NA20281.bam')->absolute], %$chrom });
}
is_deeply \@results, \@expected, 'got correct results for fofn_with_genome_chunking group_all method with chunk_size 0';

# vrpipe genome chunking with all the methods
{
    # create a fofn_with_metadata ds setup first
    ok my $ds = VRPipe::DataSource->create(
        type   => 'fofn_with_metadata',
        method => 'all',
        source => file(qw(t data calling_datasource.fofn))->absolute->stringify
      ),
      'could create a fofn_with_metadata datasource for use in next test';
    
    my $single_step = VRPipe::Step->create(
        name              => 'bam_symlink',
        inputs_definition => { bams => VRPipe::StepIODefinition->create(type => 'bam', description => 'bams', max_files => -1) },
        body_sub          => sub {
            my $self = shift;
            foreach my $bam (@{ $self->inputs->{bams} || [] }) {
                my $ofile = $self->output_file(output_key => 'the_output', basename => $bam->basename, type => 'bam');
                $bam->symlink($ofile);
            }
            return 1;
        },
        outputs_definition => { the_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'the output') },
        post_process_sub   => sub          { return 1; },
        description        => 'a step'
    );
    my $single_step_pipeline = VRPipe::Pipeline->create(name => 'bam symlink pipeline', description => 'symlink bams');
    $single_step_pipeline->add_step($single_step);
    VRPipe::StepAdaptor->create(pipeline => $single_step_pipeline, to_step => 1, adaptor_hash => { bams => { data_element => 0 } });
    my $output_root = get_output_dir('bam_symlink_output');
    my $ps = VRPipe::PipelineSetup->create(name => 'bam symlink ps', datasource => $ds, output_root => $output_root, pipeline => $single_step_pipeline, active => 0);
    $ps->trigger();
    
    my @parent_element_ids;
    my @expected_bams;
    my %symlink_path_to_parent_element_id;
    foreach my $element (@{ get_elements($ds) }) {
        my $parent_id = $element->id;
        push(@parent_element_ids, $parent_id);
        
        my ($orig_file)  = @{ $element->files() };
        my ($symlink)    = VRPipe::File->search({ parent => $orig_file->id });
        my $symlink_path = $symlink->path->stringify;
        push(@expected_bams, $symlink_path);
        $symlink_path_to_parent_element_id{$symlink_path} = $parent_id;
    }
    @parent_element_ids = sort { $a <=> $b } @parent_element_ids;
    
    # group_by_metadata
    ok $ds = VRPipe::DataSource->create(
        type    => 'vrpipe_with_genome_chunking',
        method  => 'group_by_metadata',
        source  => $ps->id . '[1]',
        options => { metadata_keys => 'sample', %$chunking_ds_options }
      ),
      'could create a vrpipe_with_genome_chunking datasource with group_by_metadata method';
    
    @results = ();
    my $correct_parents = 0;
    
    my $check_results = sub {
        my $ds = shift;
        foreach my $element (@{ get_elements($ds) }) {
            push(@results, result_with_inflated_paths($element));
            my ($symlink_path) = $element->paths();
            my $expected_parent = $symlink_path_to_parent_element_id{$symlink_path};
            my @actual_parent = VRPipe::DataElementLink->get_column_values(['parent'], { child => $element->id });
            if (@actual_parent == 1 && $actual_parent[0] == $expected_parent) {
                $correct_parents++;
            }
            else {
                $correct_parents--;
            }
        }
    };
    &$check_results($ds);
    @expected = ();
    foreach my $bam (@expected_bams) {
        my ($group) = $bam =~ /(NA\d+)\.bam$/;
        foreach my $chunk (@$chunks) {
            push(@expected, { paths => [$bam], %$chunk, group => $group });
        }
    }
    is_deeply \@results, \@expected, 'got correct results for vrpipe_with_genome_chunking group_by_metadata method';
    is $correct_parents, 42, 'all the created dataelements had the correct parent dataelements linked';
    
    # all
    ok $ds = VRPipe::DataSource->create(
        type    => 'vrpipe_with_genome_chunking',
        method  => 'all',
        source  => $ps->id . '[1]',
        options => $chunking_ds_options
      ),
      'could create a vrpipe_with_genome_chunking datasource with all method';
    
    @results         = ();
    $correct_parents = 0;
    &$check_results($ds);
    foreach my $hash (@expected) {
        delete $hash->{group};
    }
    is_deeply \@results, \@expected, 'got correct results for vrpipe_with_genome_chunking all method';
    is $correct_parents, 42, 'all the created dataelements had the correct parent dataelements linked';
    
    # group_all
    ok $ds = VRPipe::DataSource->create(
        type    => 'vrpipe_with_genome_chunking',
        method  => 'group_all',
        source  => $ps->id . '[1]',
        options => $chunking_ds_options
      ),
      'could create a vrpipe_with_genome_chunking datasource with group_all method';
    
    @results         = ();
    $correct_parents = 0;
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, result_with_inflated_paths($element));
        my $child_id = $element->id;
        my @parent_ids = sort { $a <=> $b } VRPipe::DataElementLink->get_column_values(['parent'], { child => $child_id });
        if (@parent_ids == 2 && $parent_ids[0] == $parent_element_ids[0] && $parent_ids[1] == $parent_element_ids[1]) {
            $correct_parents++;
        }
        else {
            $correct_parents--;
        }
    }
    @expected = ();
    foreach my $chunk (@$chunks) {
        push(@expected, { paths => [@expected_bams], %$chunk });
    }
    is_deeply \@results, \@expected, 'got correct results for vrpipe_with_genome_chunking group_all method';
    is $correct_parents, 21, 'all the created dataelements had the correct parent dataelements linked';
}

# confirm that a step that accepts multiple file types from the datasource only
# receives the bam files for the bam input and nothing for the other
{
    my $single_step = VRPipe::Step->create(
        name              => '2_type_step',
        inputs_definition => {
            bams => VRPipe::StepIODefinition->create(type => 'bam', description => 'bams', max_files => -1),
            vcf  => VRPipe::StepIODefinition->create(type => 'txt', min_files   => 0,      max_files => 1, description => 'a vcf')
        },
        body_sub => sub {
            my $self  = shift;
            my @bams  = @{ $self->inputs->{bams} || [] };
            my @vcfs  = @{ $self->inputs->{vcf} || [] };
            my $ofile = $self->output_file(output_key => 'the_output', basename => 'out.o', type => 'txt', metadata => { bams => scalar(@bams), vcfs => scalar(@vcfs) });
            my $fh    = $ofile->openw;
            print $fh "foo\n";
            $ofile->close;
            return 1;
        },
        outputs_definition => { the_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'the output', metadata => { bams => 'num bams', vcfs => 'num vcfs' }) },
        post_process_sub   => sub          { return 1; },
        description        => '2_type_step'
    );
    
    my $single_step_pipeline = VRPipe::Pipeline->create(name => '2_type_step pipeline', description => '2_type_step pipeline');
    $single_step_pipeline->add_step($single_step);
    VRPipe::StepAdaptor->create(pipeline => $single_step_pipeline, to_step => 1, adaptor_hash => { bams => { data_element => 0 }, vcf => { data_element => 0 } });
    
    my $output_root = get_output_dir('datasource_2type_output');
    my $ps = VRPipe::PipelineSetup->create(name => '2type_ps', datasource => $ds, output_root => $output_root, pipeline => $single_step_pipeline, active => 0);
    
    #$ps->verbose(1); # uncomment to see why this isn't working if the following test fails or the test script stalls here
    $ps->trigger();
    
    my @ss = VRPipe::StepState->search({ pipelinesetup => $ps->id, complete => 1 });
    is scalar(@ss), 21, 'a step that takes 2 file types completed when given multiple of the required type and 0 of the optional';
    my ($bams, $vcfs) = (0, 0);
    foreach my $ss (@ss) {
        my ($file) = $ss->output_files_list;
        my $meta = $file->metadata;
        $bams += $meta->{bams};
        $vcfs += $meta->{vcfs};
    }
    is_deeply [$bams, $vcfs], [63, 0], 'the step saw the correct number of input bams and vcfs';
    
    # also confirm it works when we also supply a vcf file
    $ds = VRPipe::DataSource->create(
        type    => 'fofn_with_genome_chunking',
        method  => 'group_all',
        source  => file(qw(t data bams_and_vcf.fofn))->absolute->stringify,
        options => { reference_index => $fai, chunk_override_file => $override, chrom_list => '11 20', chunk_size => 10000000 }
    );
    
    $ps = VRPipe::PipelineSetup->create(name => '2type_ps2', datasource => $ds, output_root => $output_root, pipeline => $single_step_pipeline, active => 0);
    
    #$ps->verbose(1);
    $ps->trigger();
    
    @ss = VRPipe::StepState->search({ pipelinesetup => $ps->id, complete => 1 });
    is scalar(@ss), 21, 'a step that takes 2 file types completed when given multiple of the required type and 1 of the optional';
    ($bams, $vcfs) = (0, 0);
    foreach my $ss (@ss) {
        my ($file) = $ss->output_files_list;
        my $meta = $file->metadata;
        $bams += $meta->{bams};
        $vcfs += $meta->{vcfs};
    }
    is_deeply [$bams, $vcfs], [63, 21], 'the step saw the correct number of input bams and vcfs';
}

# genome chunking with ploidy
# {
#     default=>2,
#     X =>
#     [
#         # The pseudoautosomal regions 60,001-2,699,520 and 154,931,044-155,270,560 with the ploidy 2
#         { from=>1, to=>60_000, M=>1 },
#         { from=>2_699_521, to=>154_931_043, M=>1 },
#     ],
#     Y =>
#     [
#         # No chrY in females and one copy in males
#         { from=>1, to=>59_373_566, M=>1, F=>0 },
#     ],
#     MT =>
#     [
#         # Haploid MT in males and females
#         { from=>1, to => 16_569, M=>1, F=>1 },
#     ],
#}

$chunks = [{ chrom => 'X', from => 1, to => 60000, male_ploidy => 1, female_ploidy => 2, seq_no => 1, chunk_override_file => $override }, { chrom => 'X', from => 60001, to => 2699520, male_ploidy => 2, female_ploidy => 2, seq_no => 2, chunk_override_file => $override }, { chrom => 'X', from => 2699521, to => 154931043, male_ploidy => 1, female_ploidy => 2, seq_no => 3, chunk_override_file => $override }, { chrom => 'X', from => 154931044, to => 155270560, male_ploidy => 2, female_ploidy => 2, seq_no => 4, chunk_override_file => $override }, { chrom => 'Y', from => 1, to => 59373566, male_ploidy => 1, female_ploidy => 0, seq_no => 5, chunk_override_file => $override }, { chrom => 'MT', from => 1, to => 16569, male_ploidy => 1, female_ploidy => 1, seq_no => 6, chunk_override_file => $override }];

ok $ds = VRPipe::DataSource->create(
    type    => 'fofn_with_genome_chunking',
    method  => 'all',
    source  => file(qw(t data datasource.fofn))->absolute->stringify,
    options => { reference_index => $fai, chunk_override_file => $override, chrom_list => 'X Y MT', chunk_size => 155000000, ploidy => file(qw(t data ploidy_definition))->absolute->stringify }
  ),
  'could create a fofn_with_genome_chunking datasource with all method';

@results = ();
foreach my $element (@{ get_elements($ds) }) {
    push(@results, result_with_inflated_paths($element));
}
@expected = ();
foreach my $file (file('t', 'data', 'file.bam')->absolute, file('t', 'data', 'file.cat')->absolute, file('t', 'data', 'file.txt')->absolute) {
    foreach my $chunk (@$chunks) {
        push @expected, { paths => [$file], %$chunk };
    }
}
is_deeply \@results, \@expected, 'got correct results for fofn_with_genome_chunking all method when using the ploidy option';

# test that chained vrpipe datasources update properly as they complete
{
    # first create a fofn ds and setup that a vrpipe ds can be based on
    my $single_step = VRPipe::Step->create(
        name              => 'a_step',
        inputs_definition => { the_input => VRPipe::StepIODefinition->create(type => 'any', description => 'the input', max_files => 2) },
        body_sub          => sub {
            my $self = shift;
            my ($input, $input2) = @{ $self->inputs->{the_input} };
            my $basename = $input->basename;
            my ($type) = $basename =~ /\.(bam|cat|txt)/;
            if ($input2) {
                $basename .= '.' . $input2->basename;
            }
            my $ofile = $self->output_file(output_key => 'the_output', basename => $basename . '.o', type => 'txt', metadata => { type => $type });
            # later on, we don't want the second element of the second setup to
            # instantly complete, but everything else should
            if ($basename eq 'file.cat.o') {
                $self->dispatch(["echo foo > " . $ofile->path, VRPipe::Requirements->create(memory => 10, time => 10)]);
                return 0;
            }
            else {
                my $fh = $ofile->openw;
                print $fh "foo\n";
                $ofile->close;
                return 1;
            }
        },
        outputs_definition => { the_output => VRPipe::StepIODefinition->create(type => 'txt', description => 'the output', metadata => { type => 'the type' }) },
        post_process_sub   => sub          { return 1; },
        description        => 'a step'
    );
    
    my $single_step_pipeline = VRPipe::Pipeline->create(name => 'a_pipeline', description => 'a pipeline');
    $single_step_pipeline->add_step($single_step);
    VRPipe::StepAdaptor->create(pipeline => $single_step_pipeline, to_step => 1, adaptor_hash => { the_input => { data_element => 0 } });
    
    my $fofn_ds = VRPipe::DataSource->create(
        type    => 'fofn',
        method  => 'all',
        source  => file(qw(t data datasource.fofn))->absolute->stringify,
        options => {}
    );
    
    my $output_root = get_output_dir('datasource_ps1_output');
    my $ps1 = VRPipe::PipelineSetup->create(name => 'ps1', datasource => $fofn_ds, output_root => $output_root, pipeline => $single_step_pipeline, active => 0);
    
    get_elements($fofn_ds);
    is $fofn_ds->_changed_marker, 'a7b73b4704ae4e75ebd94cc9ab43141a', 'fofn changed marker got set after elements call';
    
    # now make a vrpipe ds and setup based on the fofn one
    my $vrpipe_ds = VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'ps1[1]',
        options => {}
    );
    
    my $output_root2 = get_output_dir('datasource_ps2_output');
    my $ps2 = VRPipe::PipelineSetup->create(name => 'ps2', datasource => $vrpipe_ds, output_root => $output_root2, pipeline => $single_step_pipeline, active => 0);
    
    is $vrpipe_ds->_changed_marker, undef, 'vrpipe changed marker starts out undefined';
    my $ps2_element_count = @{ get_elements($vrpipe_ds) };
    is $vrpipe_ds->_changed_marker, '16e80f05c8d5aabea6f51165cb884958', 'vrpipe changed marker got set after elements call';
    is $ps2_element_count, 0, 'since ps1 has completed no elements, ps2 has no elements';
    
    # now make a vrpipe ds and setup that relies on both ps1 and ps2
    my $group_ds = VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'ps1[1]|ps2[1]',
        options => { metadata_keys => 'type' }
    );
    
    my $output_root3 = get_output_dir('datasource_ps3_output');
    my $ps3 = VRPipe::PipelineSetup->create(name => 'ps3', datasource => $group_ds, output_root => $output_root3, pipeline => $single_step_pipeline, active => 0);
    
    is $group_ds->_changed_marker, undef, 'group changed marker starts out undefined';
    my $ps3_element_count = @{ get_elements($group_ds) };
    is $group_ds->_changed_marker, '94bc93ad88f62c61a16b2ecc44e6702a', 'group changed marker got set after elements call';
    is $ps3_element_count, 0, 'since ps1 and ps2 have completed no elements, ps3 has no elements';
    
    # complete a single ps1 dataelement
    my ($ps1_de1, $ps1_de2, $ps1_de3) = VRPipe::DataElement->search({ datasource => $fofn_ds->id });
    $ps1->trigger(dataelement => $ps1_de1);
    get_elements($fofn_ds);
    is $fofn_ds->_changed_marker, 'a7b73b4704ae4e75ebd94cc9ab43141a', 'fofn changed marker unchanged after completing 1 ps1 element';
    $ps2_element_count = @{ get_elements($vrpipe_ds) };
    isnt $vrpipe_ds->_changed_marker, 'a7b73b4704ae4e75ebd94cc9ab43141a', 'vrpipe changed marker changed after completing 1 ps1 element';
    is $ps2_element_count, 1, 'since ps1 has completed 1 element, ps2 has 1 element';
    $ps3_element_count = @{ get_elements($group_ds) };
    my $group_change_marker = $group_ds->_changed_marker;
    isnt $group_change_marker, '94bc93ad88f62c61a16b2ecc44e6702a', 'group changed marker changed after completing 1 ps1 element';
    is $ps3_element_count, 0, 'since ps1 and ps2 are incomplete, ps3 has no elements';
    
    # complete a single ps2 dataelement
    my ($ps2_de1) = VRPipe::DataElement->search({ datasource => $vrpipe_ds->id });
    $ps2->trigger(dataelement => $ps2_de1);
    $ps2_element_count = @{ get_elements($vrpipe_ds) };
    is $vrpipe_ds->_changed_marker, 'c7329a97a6906d3fa562aa07d62fcb04', 'vrpipe changed marker unchanged after completing its element';
    is $ps2_element_count, 1, 'since ps1 has completed 1 element, ps2 has 1 element';
    $ps3_element_count = @{ get_elements($group_ds) };
    isnt $group_ds->_changed_marker, $group_change_marker, 'group changed marker changed after completing 1 ps2 element';
    is $ps3_element_count, 0, 'since ps1 and ps2 are still incomplete, ps3 has no elements';
    
    # complete all of ps1, then trigger all 3 setups simultaneously
    $ps1->trigger;
    my $fm = Parallel::ForkManager->new(3);
    foreach my $setup ($ps1, $ps2, $ps3) {
        $fm->start and next;
        $setup->trigger();
        $fm->finish(0);
    }
    $fm->wait_all_children;
    
    my @ps1_file_paths;
    $ps2_element_count = 0;
    foreach my $element (@{ get_elements($vrpipe_ds) }) {
        $ps2_element_count++;
        push(@ps1_file_paths, map { @{ $_->{paths} } } result_with_inflated_paths($element));
    }
    
    is $vrpipe_ds->_changed_marker, 'b3ed212622b776f82d266ea1d239db08', 'vrpipe changed marker changed after ps1 completed';
    is $ps2_element_count, 3, 'since ps1 has completed 3 elements, ps2 has 3 elements';
    
    # there was a bug that meant the following test would fail due to elements
    # getting created based on just ps1's files, since ps2 got ignored. Due to
    # timing issues with the parallel forks, it didn't alwas fail though; ps3
    # needs to trigger fractions of a second after ps2
    $ps3_element_count = VRPipe::DataElement->search({ datasource => $group_ds->id });
    is $ps3_element_count, 0, 'since ps2 is still incomplete, ps3 has no elements';
    
    # test the include_in_all_elements option of the group_by_metadata method
    my $iiae_ds = VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'group_by_metadata',
        source  => 'ps1[1]',
        options => { metadata_keys => 'type', include_in_all_elements => 'ps2[1]' }
    );
    
    my $output_root4 = get_output_dir('datasource_ps4_output');
    my $ps4 = VRPipe::PipelineSetup->create(name => 'ps4', datasource => $iiae_ds, output_root => $output_root4, pipeline => $single_step_pipeline, active => 0);
    
    my $ps4_element_count = @{ get_elements($iiae_ds) };
    is $ps4_element_count, 0, 'since ps2 is still incomplete, ps4 has no elements';
    my $de_to_withdraw;
    {
        # force complete setup 2
        foreach my $des (VRPipe::DataElementState->search({ pipelinesetup => $ps2->id })) {
            $des->completed_steps(1);
            $des->update;
            $de_to_withdraw ||= $des->dataelement;
        }
        foreach my $ss (VRPipe::StepState->search({ pipelinesetup => $ps2->id })) {
            $ss->complete(1);
            $ss->update;
        }
    }
    my @ps4_elements = @{ get_elements($iiae_ds) };
    is scalar(@ps4_elements), 3, 'now ps2 is complete, ps4 has expected elements';
    my (@expected_all, @expected_ps2_de_ids);
    foreach my $ss ($ps2->states) {
        push(@expected_all, map { $_->path->stringify } $ss->output_files_list);
        push(@expected_ps2_de_ids, $ss->dataelement->id);
    }
    my (@expected_paths_per_element, @expected_parents_per_element);
    foreach my $ss ($ps1->states) {
        my ($path) = map { $_->path->stringify } $ss->output_files_list;
        push(@expected_paths_per_element, [$path, @expected_all]);
        push(@expected_parents_per_element, [$ss->dataelement->id, @expected_ps2_de_ids]);
    }
    my (@actual_paths_per_element, @actual_parents_per_element);
    foreach my $element (@ps4_elements) {
        push(@actual_paths_per_element, [map { $_->path->stringify } $element->filelist->files]);
        my @parents;
        foreach my $del (VRPipe::DataElementLink->search({ child => $element->id })) {
            push(@parents, $del->parent->id);
        }
        push(@actual_parents_per_element, [sort { $a <=> $b } @parents]);
    }
    is_deeply \@actual_paths_per_element,   \@expected_paths_per_element,   'each ps4 element has all ps2 files with a single ps1 file due to include_in_all_elements';
    is_deeply \@actual_parents_per_element, \@expected_parents_per_element, 'each ps4 element has links to every element that made the input files, including the include_in_all_elements files';
    
    # test the include_in_all_elements option of the all method
    $iiae_ds = VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'ps1[1]',
        options => { include_in_all_elements => 'ps2[1]' }
    );
    
    my $output_root5 = get_output_dir('datasource_ps5_output');
    my $ps5 = VRPipe::PipelineSetup->create(name => 'ps5', datasource => $iiae_ds, output_root => $output_root5, pipeline => $single_step_pipeline, active => 0);
    
    my @ps5_elements = @{ get_elements($iiae_ds) };
    is scalar(@ps5_elements), 3, 'ps5 (all with include_in_all_elements) has expected number of elements';
    (@actual_paths_per_element, @actual_parents_per_element) = ();
    foreach my $element (@ps5_elements) {
        push(@actual_paths_per_element, [map { $_->path->stringify } $element->filelist->files]);
        my @parents;
        foreach my $del (VRPipe::DataElementLink->search({ child => $element->id })) {
            push(@parents, $del->parent->id);
        }
        push(@actual_parents_per_element, [sort { $a <=> $b } @parents]);
    }
    is_deeply \@actual_paths_per_element,   \@expected_paths_per_element,   'each ps5 element has all ps2 files with a single ps1 file due to include_in_all_elements';
    is_deeply \@actual_parents_per_element, \@expected_parents_per_element, 'each ps5 element has links to every element that made the input files, including the include_in_all_elements files';
    
    # test that the include_in_all_elements datasource update when the include
    # setup(s) change
    my $iiae_ds_changed_marker = $iiae_ds->_changed_marker;
    ok $iiae_ds_changed_marker, 'ps5 changed marker starts with some value';
    is $iiae_ds->_source_instance->_has_changed, 1, '_has_changed returns 1 due to the transition from 0 elements to non-0';
    is $iiae_ds->_source_instance->_has_changed, 0, '_has_changed returns 0 following no further changes';
    $de_to_withdraw->withdrawn(1);
    $de_to_withdraw->update;
    is $iiae_ds->_source_instance->_has_changed, 1, '_has_changed returns 1 after withdrawing one of the include_in_all_elements setup elements';
    
    # pretend ps1's datasource uses DataSourceRole's
    # _start_over_elements_due_to_file_metadata_change() and 2 files changed
    # to test that we only start_over ps5's elements once... I don't know how
    # to actually test for that without looking at debug output, hence these
    # tests are commented out but behaved as desired after solving this problem
    # $ps5->trigger();
    # my $anti_repeat_store = {};;
    # $fofn_ds->_source_instance->_start_over_elements_due_to_file_metadata_change({ changed => [[VRPipe::File->get(path => file(qw(t data file.bam))->absolute)]] }, ['foo', 'bar'], $anti_repeat_store);
    # $fofn_ds->_source_instance->_start_over_elements_due_to_file_metadata_change({ changed => [[VRPipe::File->get(path => file(qw(t data file.cat))->absolute)]] }, ['foo', 'bar'], $anti_repeat_store);
    
    # test the filter and graph_filter options of the vrpipe datasource
    my $schema    = VRPipe::Schema->create("VRTrack");
    my @ps1_nodes = map { $schema->add_file($_) } @ps1_file_paths;
    my $vrsample1 = $schema->add("Sample", { name => "sample1", qc_failed => 0 }, outgoing => { node => $ps1_nodes[0], type => 'has' });
    my $vrsample2 = $schema->add("Sample", { name => "sample2", qc_failed => 1 }, outgoing => { node => $ps1_nodes[1], type => 'has' });
    my $vrsample3 = $schema->add("Sample", { name => "sample3" }, outgoing => { node => $ps1_nodes[2], type => 'has' });
    my @ps1_files = map { VRPipe::File->get(path => $_) } @ps1_file_paths;
    $ps1_files[1]->add_metadata({ filtkey => 'filtvalue' });
    my $filt_ds = VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => 'ps1[1]',
        options => { filter => 'filtkey#filtvalue', graph_filter => 'VRTrack#Sample#qc_failed#0' }
    );
    my @filt_elements = @{ get_elements($filt_ds) };
    is scalar(@filt_elements), 0, 'filter and graph_filter options can cancel each other out';
    $ps1_files[2]->add_metadata({ filtkey => 'filtvalue' });
    @filt_elements = @{ get_elements($filt_ds) };
    is scalar(@filt_elements), 1, 'filter and graph_filter options can work together correctly';
    $ps1_files[0]->add_metadata({ filtkey => 'filtvalue' });
    @filt_elements = @{ get_elements($filt_ds) };
    is scalar(@filt_elements), 2, 'graph_filter with a 0 value can get both 0 value and unset properties';
}

# author-only tests for the irods and graph_vrtrack datasources
SKIP: {
    my $num_tests = 60;
    skip "author-only tests for an iRods datasource", $num_tests unless ($ENV{VRPIPE_AUTHOR_TESTS} && $ENV{VRPIPE_IRODS_TEST_ROOT} && $ENV{VRPIPE_IRODS_TEST_RESOURCE});
    
    my $output_root    = get_output_dir('datasource_irods_import_dir');
    my $irods_root     = $ENV{VRPIPE_IRODS_TEST_ROOT};
    my $irods_resource = $ENV{VRPIPE_IRODS_TEST_RESOURCE};
    my (undef, $irods_zone) = split('/', $irods_root);
    my $schema = VRPipe::Schema->create("VRTrack");
    
    # irods ds uses IPC::Run, which breaks if STDERR is tied, so we'll tie it
    # to test that our workaround works
    {
        package MySTDERR;
        use Symbol qw(geniosym);
        sub TIEHANDLE { return bless geniosym, __PACKAGE__ }
        sub PRINT { shift; print @_ }
    }
    tie *STDERR, 'MySTDERR' or die $!;
    
    # real-world test
    ok my $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'all',
        source  => 'seq',
        options => {
            file_query     => q[sequenom_plate LIKE '%' and study_id = 2622 and dcterms:created '<' 2013-07-26],
            local_root_dir => $output_root
        }
      ),
      'could create an irods datasource';
    
    my @results;
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, result_with_inflated_paths($element));
    }
    is_deeply \@results, [{ paths => [file($output_root, qw(seq sequenom 05 94 43 QC288261____20130701_G01.csv))], irods_path => '/seq/sequenom/05/94/43/QC288261____20130701_G01.csv' }, { paths => [file($output_root, qw(seq sequenom 14 62 84 QC288261____20130701_C01.csv))], irods_path => '/seq/sequenom/14/62/84/QC288261____20130701_C01.csv' }, { paths => [file($output_root, qw(seq sequenom 95 35 0e QC288261____20130701_A01.csv))], irods_path => '/seq/sequenom/95/35/0e/QC288261____20130701_A01.csv' }, { paths => [file($output_root, qw(seq sequenom d8 7c 21 QC288261____20130701_E01.csv))], irods_path => '/seq/sequenom/d8/7c/21/QC288261____20130701_E01.csv' }], 'got correct results for irods all';
    
    my $file = VRPipe::File->create(path => file($output_root, qw(seq sequenom 05 94 43 QC288261____20130701_G01.csv)));
    my $expected_file_meta = { sample_cohort => '20f8a331-69ac-4510-94ab-e3a69c50e46f', sequenom_well => 'G01', sample_common_name => 'Homo sapiens', sequenom_plate => 'QC288261____20130701', study_id => 2622, sample_consent => 1, sample_supplier_name => 'd2b57a6a-9dd8-4e7d-868e-9209a399711b', sample_id => 1653292, sample => 'QC1Hip-1', sample_accession_number => 'SAMEA2398742', sample_control => 0, md5 => '059443dbff29215ff8b6aa6e247b072f', irods_path => '/seq/sequenom/05/94/43/QC288261____20130701_G01.csv', manual_qc => 1, sequenom_plex => 'W30467', expected_md5 => '059443dbff29215ff8b6aa6e247b072f', sample_donor_id => '20f8a331-69ac-4510-94ab-e3a69c50e46f' };
    is_deeply $file->metadata, $expected_file_meta, 'correct file metadata was present on one of the irods files';
    
    # check the warehouse method
    ok $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'all_with_warehouse_metadata',
        source  => 'seq',
        options => {
            file_query        => q[sequenom_plate LIKE '%' and study_id = 2622 and dcterms:created '<' 2013-07-26],
            local_root_dir    => $output_root,
            required_metadata => 'sample_cohort,public_name'
        }
      ),
      'could create an irods datasource with all_with_warehouse_metadata method';
    is scalar(@{ get_elements($ds) }), scalar(@results), 'required_metadata option gave us the correct number of elements for the all_with_warehouse_metadata method';
    $file->reselect_values_from_db;
    $expected_file_meta->{public_name}         = 'HPSI0813i-ffdb_3';
    $expected_file_meta->{sample_created_date} = '2013-06-25 14:09:22';
    $expected_file_meta->{taxon_id}            = 9606;
    $expected_file_meta->{study_title}         = 'G0325 [collection qc1] Wellcome Trust Strategic Award application – HIPS';
    is_deeply $file->metadata, $expected_file_meta, 'correct file metadata was present on one of the irods files, including warehouse metadata';
    
    my $output_root_grouped = get_output_dir('datasource_irods_import_dir_grouped');
    ok $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'group_by_metadata_with_warehouse_metadata',
        source  => 'seq',
        options => {
            file_query        => q[sequenom_plate LIKE '%' and study_id = 2622 and dcterms:created '<' 2013-07-26],
            local_root_dir    => $output_root_grouped,
            required_metadata => 'sample_cohort,public_name',
            metadata_keys     => 'study_id'
        }
      ),
      'could create an irods datasource with group_by_metadata_with_warehouse_metadata method';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, result_with_inflated_paths($element));
    }
    is_deeply \@results, [{ paths => [file($output_root_grouped, qw(seq sequenom 05 94 43 QC288261____20130701_G01.csv)), file($output_root_grouped, qw(seq sequenom 14 62 84 QC288261____20130701_C01.csv)), file($output_root_grouped, qw(seq sequenom 95 35 0e QC288261____20130701_A01.csv)), file($output_root_grouped, qw(seq sequenom d8 7c 21 QC288261____20130701_E01.csv))], group => '2622' }], 'got correct results for irods group_by_metadata_with_warehouse_metadata';
    is_deeply [VRPipe::File->get(path => file($output_root_grouped, qw(seq sequenom 05 94 43 QC288261____20130701_G01.csv)))->metadata->{irods_path}, VRPipe::File->get(path => file($output_root_grouped, qw(seq sequenom 14 62 84 QC288261____20130701_C01.csv)))->metadata->{irods_path}, VRPipe::File->get(path => file($output_root_grouped, qw(seq sequenom 95 35 0e QC288261____20130701_A01.csv)))->metadata->{irods_path}, VRPipe::File->get(path => file($output_root_grouped, qw(seq sequenom d8 7c 21 QC288261____20130701_E01.csv)))->metadata->{irods_path}], ['/seq/sequenom/05/94/43/QC288261____20130701_G01.csv', '/seq/sequenom/14/62/84/QC288261____20130701_C01.csv', '/seq/sequenom/95/35/0e/QC288261____20130701_A01.csv', '/seq/sequenom/d8/7c/21/QC288261____20130701_E01.csv'], 'files have correct irods_path metadata for group_by_metadata_with_warehouse_metadata';
    
    # check that group_by_metadata_with_warehouse_metadata works with no
    # local_root_dir
    ok $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'group_by_metadata_with_warehouse_metadata',
        source  => 'seq',
        options => {
            file_query    => q[target = 1 and study_id = 3474 and manual_qc = 1 and type = cram and sample = vbseqx106041591],
            metadata_keys => 'library'
        }
      ),
      'could create an irods datasource with group_by_metadata_with_warehouse_metadata method and no local_root_dir';
    
    @results = ();
    my $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/16139/16139_7.cram', 'irods:/seq/16139/16139_8.cram'], 1], 'the resulting files had irods protocols, and grouping working correctly';
    
    # also check that it works with STDERR untied
    untie *STDERR;
    
    # check that we can aggregate results from multiple imeta queries specified
    # in a file
    ok $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'all_with_warehouse_metadata',
        source  => 'archive',
        options => {
            file_query     => file(qw(t data irods.datasource))->absolute->stringify,
            local_root_dir => $output_root
        }
      ),
      'could create an irods datasource with all_with_warehouse_metadata method';
    my $els = get_elements($ds);
    is scalar(@$els), 27, 'aggregating 2 imeta queries using a file source worked';
    
    # check the warehouse method with the option to get qc-related files, and
    # also check we don't need a local_root_dir
    ok $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'all_with_warehouse_metadata',
        source  => 'seq',
        options => {
            file_query       => q[id_run = 15744 and target = 1 and type = cram],
            require_qc_files => 1
        }
      ),
      'could create an irods datasource using all_with_warehouse_metadata method with require_qc_files option and no local_root_dir';
    is scalar(@{ get_elements($ds) }), 8, 'got the correct number of elements';
    $file = VRPipe::File->get(path => '/seq/15744/15744_8.cram', protocol => 'irods:');
    $expected_file_meta = { expected_md5 => 'af0d56cce970925d7bba207784a5e1c2', md5 => 'af0d56cce970925d7bba207784a5e1c2', reference => '/lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa0_6/hs37d5.fa', lane => '15744_8', library => 'HiSeqX_NX_Titration_NA19239_H 13237756', study => 'HX Test Plan', id_run => 15744, library_id => 13237756, sample_id => 2247346, target => 1, study_title => 'HX Test Plan', total_reads => 623053650, alignment => 1, study_id => 3165, sample => 'HiSeqX_NX_Titration_NA19239_H', sample_created_date => '2015-03-12 08:49:13', is_paired_read => 1 };
    is_deeply $file->metadata, $expected_file_meta, 'correct file metadata was present on one of the irods files';
    ok my $graph_file = $schema->get_file($file->protocolless_path->stringify, $file->protocol), 'there was a node in the graph db for one of the cram files';
    my @qc_files = $graph_file->related(outgoing => { type => 'qc_file' }) if $graph_file;
    is_deeply [map { $_->path } sort { $a->path cmp $b->path } @qc_files], ['irods:/seq/15744/15744_8_F0xB00.stats', 'irods:/seq/15744/qc/15744_8.genotype.json', 'irods:/seq/15744/qc/15744_8.verify_bam_id.json'], 'irods qc files were associated with the cram file';
    
    # more complete test with our own freshly-added files and metadata
    system("irm -fr $irods_root > /dev/null 2> /dev/null");
    system("imkdir -p $irods_root");
    system("iput -R $irods_resource t/data/file.txt $irods_root");
    system("iput -R $irods_resource t/data/file2.txt $irods_root");
    system("imeta -z $irods_zone add -d $irods_root/file.txt study_id 2623");
    system("imeta -z $irods_zone add -d $irods_root/file2.txt study_id 2623");
    system("imeta -z $irods_zone add -d $irods_root/file.txt foo bar");
    system("imeta -z $irods_zone add -d $irods_root/file2.txt simon says");
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'all',
        source  => $irods_zone,
        options => {
            file_query      => q[study_id = 2623],
            local_root_dir  => $output_root,
            update_interval => 10
        }
      ),
      'could create another irods datasource';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, result_with_inflated_paths($element));
    }
    
    is_deeply \@results, [{ paths => [file($output_root, $irods_root, 'file.txt')], irods_path => "$irods_root/file.txt" }, { paths => [file($output_root, $irods_root, 'file2.txt')], irods_path => "$irods_root/file2.txt" }], 'got correct results for irods all';
    
    $file = VRPipe::File->create(path => file($output_root, $irods_root, 'file.txt'));
    is_deeply $file->metadata, { study_id => 2623, foo => 'bar', irods_path => "$irods_root/file.txt" }, 'correct file metadata was present on one of the irods files';
    
    # alter a bit of metadata and test that VRPipe notices the change
    system("imeta -z $irods_zone mod -d $irods_root/file.txt foo bar v:car");
    get_elements($ds);
    $file->reselect_values_from_db;
    is_deeply $file->metadata, { study_id => 2623, foo => 'bar', irods_path => "$irods_root/file.txt" }, 'metadata in VRPipe unchanged after change in irods, due to cached result';
    sleep(11);
    get_elements($ds);
    $file->reselect_values_from_db;
    is_deeply $file->metadata, { study_id => 2623, foo => 'car', irods_path => "$irods_root/file.txt" }, 'metadata in VRPipe updated correctly after waiting 10 seconds';
    
    # add a new file to make sure we pick that up as well
    system("iput -R $irods_resource t/data/file3.txt $irods_root");
    system("imeta -z $irods_zone add -d $irods_root/file3.txt study_id 2623");
    system("imeta -z $irods_zone add -d $irods_root/file3.txt simple simon");
    sleep(11);
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, result_with_inflated_paths($element));
    }
    
    my @local_files = (file($output_root, $irods_root, 'file.txt'), file($output_root, $irods_root, 'file2.txt'), file($output_root, $irods_root, 'file3.txt'));
    is_deeply \@results, [{ paths => [$local_files[0]], irods_path => "$irods_root/file.txt" }, { paths => [$local_files[1]], irods_path => "$irods_root/file2.txt" }, { paths => [$local_files[2]], irods_path => "$irods_root/file3.txt" }], 'got correct results for irods all after adding a new file';
    
    $file = VRPipe::File->create(path => $local_files[2]);
    is_deeply $file->metadata, { study_id => 2623, simple => 'simon', irods_path => "$irods_root/file3.txt" }, 'correct file metadata was present on the newly added irods file';
    
    # test graph_filter option of the irods datasource
    system("imeta -z $irods_zone add -d $irods_root/file.txt sample samplea");
    system("imeta -z $irods_zone add -d $irods_root/file2.txt sample sampleb");
    system("imeta -z $irods_zone add -d $irods_root/file3.txt sample sampleb");
    my $filt_ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'all_with_warehouse_metadata',
        source  => $irods_zone,
        options => {
            file_query      => q[study_id = 2623],
            local_root_dir  => $output_root,
            update_interval => 10,
            graph_filter    => 'VRTrack#Sample#qc_failed#0'
        }
    );
    my @filt_elements = @{ get_elements($filt_ds) };
    is scalar(@filt_elements), 3, 'all_with_warehouse_metadata graph_filter with a 0 value passes all files with the property unset';
    my @lf_nodes = map { $schema->add_file($_->stringify, 'irods:') } @local_files;
    my $vrsamplea = $schema->get("Sample", { name => "samplea" });
    $vrsamplea->qc_failed(0);
    my $vrsampleb = $schema->add("Sample", { name => "sampleb" });
    $vrsampleb->qc_failed(1);
    sleep(11);
    @filt_elements = @{ get_elements($filt_ds) };
    is scalar(@filt_elements), 1, 'changes to what would pass the graph_filter are picked up correctly';
    
    $filt_ds = VRPipe::DataSource->create(
        type    => 'irods',
        method  => 'group_by_metadata_with_warehouse_metadata',
        source  => $irods_zone,
        options => {
            file_query      => q[study_id = 2623],
            local_root_dir  => $output_root,
            update_interval => 10,
            graph_filter    => 'VRTrack#Sample#qc_failed#1',
            metadata_keys   => 'sample'
        }
    );
    @filt_elements = @{ get_elements($filt_ds) };
    @results = map { result_with_inflated_paths($_) } @filt_elements;
    is_deeply [scalar(@filt_elements), $results[0]], [1, { paths => [$local_files[1], $local_files[2]], group => "sampleb" }], 'group_by_metadata_with_warehouse_metadata graph_filter with a 1 value returns a single dataelement with the 2 paths associated with the passing sample';
    
    system("irm -fr $irods_root");
    
    # after running the irods datasources we have study 3165 with multiple
    # samples and cram files that we can test the graph_vrtrack datasource on,
    # but for a full test we need to fake that we've run npg_cram_parser step
    # by manually adding simplified Bam_Stats and Genotype nodes to 2 of the
    # crams
    my @fake_stats = ({ a => 99, b => 10001, c => 0, d => 1 }, { a => 101, b => 9999, c => 1, d => 0 });
    foreach my $cram_path ('/seq/15744/15744_6.cram', '/seq/15744/15744_4.cram') {
        my $cram_node = $schema->get_file($cram_path, 'irods:');
        my $stats = shift(@fake_stats);
        
        foreach my $qc_file ($cram_node->related(outgoing => { type => 'qc_file' })) {
            my $qc_path = $qc_file->path;
            
            if ($qc_path =~ /\.stats$/) {
                $schema->add(
                    'Bam_Stats',
                    {
                        mode                  => 'normal',
                        options               => '-opt',
                        date                  => 1443015733,
                        'reads QC failed'     => $stats->{a},
                        'raw total sequences' => $stats->{b}
                    },
                    incoming => { type => 'summary_stats', node => $qc_file }
                );
            }
            elsif ($qc_path =~ /\.genotype\.json$/) {
                $schema->add(
                    'Genotype',
                    {
                        date                 => 1443015733,
                        pass                 => $stats->{c},
                        expected_sample_name => 'foo',
                        matched_sample_name  => 'foo'
                    },
                    incoming => { type => 'genotype_data', node => $qc_file }
                );
            }
            elsif ($qc_path =~ /\.verify_bam_id\.json$/) {
                $schema->add(
                    'Verify_Bam_ID',
                    {
                        date    => 1443015733,
                        pass    => $stats->{d},
                        freemix => 'foo'
                    },
                    incoming => { type => 'verify_bam_id_data', node => $qc_file }
                );
            }
        }
    }
    
    # to test a parent_filter we'll also set some qc fail info on some samples
    @fake_stats = (1, 0, 1, undef, 1, 0, 1, 0);
    my $study_node = $schema->get('Study', { id => 3165 });
    foreach my $sample_node (sort { $a->id <=> $b->id } $study_node->related(outgoing => { type => 'member' })) {
        my $failed = shift(@fake_stats);
        next unless defined $failed;
        $sample_node->qc_failed($failed);
    }
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => {}
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and no other options';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_1.cram', 'irods:/seq/15744/15744_2.cram', 'irods:/seq/15744/15744_3.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_5.cram', 'irods:/seq/15744/15744_6.cram', 'irods:/seq/15744/15744_7.cram', 'irods:/seq/15744/15744_8.cram'], 8], 'the correct irods files were returned in the correct dataelements';
    
    $expected_file_meta = {
        gender_gender       => 'U',
        gender_source       => 'sequencescape',
        lane_lane           => 1,
        sample_qc_failed    => 1,
        alignment           => 1,
        expected_md5        => 'b02005ae2e76ed683741c72f174ff636',
        id_run              => 15744,
        is_paired_read      => 1,
        lane                => '15744_1',
        library             => 'HiSeqX_NX_Titration_NA19239_A 13237749',
        library_id          => 13237749,
        md5                 => 'b02005ae2e76ed683741c72f174ff636',
        reads               => 457093492,
        reference           => '/lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa0_6/hs37d5.fa',
        sample              => 'HiSeqX_NX_Titration_NA19239_A',
        sample_created_date => '2015-03-12 08:49:11',
        sample_id           => 2247339,
        study               => 'HX Test Plan',
        study_id            => 3165,
        study_title         => 'HX Test Plan',
        target              => 1,
        total_reads         => 457093492
    };
    my $example_cram_vrfile = VRPipe::File->get(path => '/seq/15744/15744_1.cram', protocol => 'irods:');
    is_deeply $example_cram_vrfile->metadata, $expected_file_meta, 'the metadata on one of the cram files was correct';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { parent_filter => 'Sample#qc_failed#0' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a parent filter';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_2.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_6.cram', 'irods:/seq/15744/15744_8.cram'], 4], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { parent_filter => 'Sample#qc_failed#0,Sample#created_date#1426150152' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and 2 parent filters';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_2.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_6.cram'], 3], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { qc_filter => 'stats#raw total sequences#>#10000' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a stats qc filter';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_6.cram'], 1], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { qc_filter => 'genotype#pass#=#1' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a genotype qc filter';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_4.cram'], 1], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { qc_filter => 'verifybamid#pass#=#1' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a verifybamid qc filter';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_6.cram'], 1], 'the correct irods files were returned in the correct dataelements';
    
    $example_cram_vrfile = $schema->get_file('/seq/15744/15744_2.cram', 'irods:');
    $example_cram_vrfile->add_properties({ target => 0 });
    $example_cram_vrfile = $schema->get_file('/seq/15744/15744_7.cram', 'irods:');
    $example_cram_vrfile->remove_property('target');
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { qc_filter => 'file#target#=#1' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a file qc filter';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_1.cram', 'irods:/seq/15744/15744_3.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_5.cram', 'irods:/seq/15744/15744_6.cram', 'irods:/seq/15744/15744_8.cram'], 6], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { qc_filter => 'stats#raw total sequences#>#10000,stats#reads QC failed#<#100,genotype#pass#=#1,verifybamid#pass#=#1,file#target#=#1' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and multiple qc filters that together match no files';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [[], 0], 'no irods files were returned in zero dataelements, as expected';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { qc_filter => 'stats#raw total sequences#>#9998,stats#reads QC failed#<#102,genotype#pass#=#1,verifybamid#pass#=#0,file#target#=#1' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and multiple qc filters that together match 1 file';
    
    @results       = ();
    $element_count = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $element_count++;
        push(@results, $element->paths);
    }
    is_deeply [\@results, $element_count], [['irods:/seq/15744/15744_4.cram'], 1], 'the correct irods files were returned in the correct dataelements';
    
    # metadata corresponding to what we filtered on should be set on the file
    $expected_file_meta = {
        'study_id'                  => '3165',
        'reads'                     => '787414252',
        'gender_source'             => 'sequencescape',
        'is_paired_read'            => '1',
        'library_id'                => '13237752',
        'library'                   => 'HiSeqX_NX_Titration_NA19239_D 13237752',
        'study_title'               => 'HX Test Plan',
        'stats_raw total sequences' => '9999',
        'target'                    => '1',
        'reference'                 => '/lustre/scratch109/srpipe/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa0_6/hs37d5.fa',
        'alignment'                 => '1',
        'gender_gender'             => 'U',
        'sample'                    => 'HiSeqX_NX_Titration_NA19239_D',
        'expected_md5'              => 'bd229470c44887e9a19f2a41533d5c8e',
        'study'                     => 'HX Test Plan',
        'lane'                      => '15744_4',
        'sample_created_date'       => '2015-03-12 08:49:12',
        'verifybamid_pass'          => '0',
        'total_reads'               => '787414252',
        'id_run'                    => '15744',
        'stats_reads QC failed'     => '101',
        'lane_lane'                 => '4',
        'genotype_pass'             => '1',
        'sample_id'                 => '2247342',
        'md5'                       => 'bd229470c44887e9a19f2a41533d5c8e'
    };
    $example_cram_vrfile = VRPipe::File->get(path => '/seq/15744/15744_4.cram', protocol => 'irods:');
    is_deeply $example_cram_vrfile->metadata, $expected_file_meta, 'the metadata on one of the cram files was correct';
    
    my $male = $schema->get('Gender', { gender => 'M' });
    foreach my $path ('/seq/15744/15744_2.cram', '/seq/15744/15744_3.cram', '/seq/15744/15744_6.cram') {
        my $cram = $schema->get_file($path, 'irods:');
        my $sample = $cram->closest('VRTrack', 'Sample', direction => 'incoming');
        $sample->relate_to($male, 'gender', replace => 1);
    }
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { group_by_metadata => 'Sample' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a group on sample';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, [$element->paths]);
    }
    is_deeply \@results, [['irods:/seq/15744/15744_1.cram'], ['irods:/seq/15744/15744_2.cram'], ['irods:/seq/15744/15744_3.cram'], ['irods:/seq/15744/15744_4.cram'], ['irods:/seq/15744/15744_5.cram'], ['irods:/seq/15744/15744_6.cram'], ['irods:/seq/15744/15744_7.cram'], ['irods:/seq/15744/15744_8.cram']], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { group_by_metadata => 'Study' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a group on study';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, [$element->paths]);
    }
    is_deeply \@results, [['irods:/seq/15744/15744_1.cram', 'irods:/seq/15744/15744_2.cram', 'irods:/seq/15744/15744_3.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_5.cram', 'irods:/seq/15744/15744_6.cram', 'irods:/seq/15744/15744_7.cram', 'irods:/seq/15744/15744_8.cram']], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { group_by_metadata => 'Gender' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a group on gender';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, [$element->paths]);
    }
    is_deeply \@results, [['irods:/seq/15744/15744_1.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_5.cram', 'irods:/seq/15744/15744_7.cram', 'irods:/seq/15744/15744_8.cram'], ['irods:/seq/15744/15744_2.cram', 'irods:/seq/15744/15744_3.cram', 'irods:/seq/15744/15744_6.cram']], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { group_by_metadata => 'Study,Gender' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a group on study and gender';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, [$element->paths]);
    }
    is_deeply \@results, [['irods:/seq/15744/15744_1.cram', 'irods:/seq/15744/15744_4.cram', 'irods:/seq/15744/15744_5.cram', 'irods:/seq/15744/15744_7.cram', 'irods:/seq/15744/15744_8.cram'], ['irods:/seq/15744/15744_2.cram', 'irods:/seq/15744/15744_3.cram', 'irods:/seq/15744/15744_6.cram']], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => { group_by_metadata => 'Sample,Gender' }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method and a group on sample and gender';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, [$element->paths]);
    }
    is_deeply \@results, [['irods:/seq/15744/15744_1.cram'], ['irods:/seq/15744/15744_2.cram'], ['irods:/seq/15744/15744_3.cram'], ['irods:/seq/15744/15744_4.cram'], ['irods:/seq/15744/15744_5.cram'], ['irods:/seq/15744/15744_6.cram'], ['irods:/seq/15744/15744_7.cram'], ['irods:/seq/15744/15744_8.cram']], 'the correct irods files were returned in the correct dataelements';
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'graph_vrtrack',
        method  => 'lanelet_crams',
        source  => 'Study#id#3165',
        options => {
            group_by_metadata => 'Gender',
            qc_filter         => 'stats#raw total sequences#>#9998,stats#reads QC failed#<#102,file#target#=#1'
        }
      ),
      'could create a graph_vrtrack datasource with lanelet_crams method, grouped on gender and with multiple filters';
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, [$element->paths]);
    }
    is_deeply \@results, [['irods:/seq/15744/15744_4.cram'], ['irods:/seq/15744/15744_6.cram']], 'the correct irods files were returned in the correct dataelements';
}

# test a special vrtrack test database; these tests are meant for the author
# only, but will also work for anyone with a working VRTrack::Factory setup
# (the '_has_changed gives no change' tests were written at a time when
#  _has_changed did not set _changed_marker, so are a little weird now but were
#  kept anyway)
SKIP: {
    my $num_tests = 31;
    skip "author-only tests for a VRTrack datasource", $num_tests unless $ENV{VRPIPE_VRTRACK_TESTDB};
    eval "require VRTrack::Factory;";
    skip "VRTrack::Factory not loading", $num_tests if $@;
    
    # create the vrtrack db
    my %cd  = VRTrack::Factory->connection_details('rw');
    my @sql = VRTrack::VRTrack->schema();
    open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
    print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    foreach my $sql (@sql) {
        print $mysqlfh $sql;
    }
    close($mysqlfh);
    
    # populate it *** uses update_vrmeta.pl, which may not be supplied with
    #                 VRTrack in the future...
    system("update_vrmeta.pl --samples t/data/vrtrack.samples --index t/data/vrtrack.sequence.index --database $ENV{VRPIPE_VRTRACK_TESTDB} > /dev/null 2> /dev/null");
    
    # alter things on the lanes to enable useful tests
    my $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $lanes = $vrtrack->processed_lane_hnames();
    my %expectations;
    for my $i (1 .. 60) {
        my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lanes->[$i - 1]);
        my $name = $lane->hierarchy_name;
        
        if ($i <= 50) {
            $lane->auto_qc_status('passed');
            push(@{ $expectations{auto_qc_passed} }, $name);
            if ($i > 10) {
                $lane->genotype_status('confirmed');
                push(@{ $expectations{gt_status_confirmed} }, $name);
            }
            if ($i > 20) {
                $lane->qc_status('pending');
                push(@{ $expectations{pending} }, $name);
            }
            if ($i > 30) {
                $lane->qc_status('no_qc');
                push(@{ $expectations{no_qc} }, $name);
            }
            if ($i > 40) {
                $lane->qc_status('failed');
                push(@{ $expectations{failed} }, $name);
            }
            if ($i > 45) {
                $lane->qc_status('passed');
                push(@{ $expectations{qc_status_passed} }, $name);
            }
        }
        elsif ($i > 55) {
            $lane->is_withdrawn(1);
            push(@{ $expectations{withdrawn} }, $name);
        }
        
        $lane->update;
    }
    
    ok $ds = VRPipe::DataSource->create(
        type    => 'vrtrack',
        method  => 'lanes',
        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
        options => { auto_qc_status => 'passed' }
      ),
      'could create a vrtrack datasource';
    my $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 50, 'got correct number of results for vrtrack lanes auto_qc_status => passed';
    
    ## tests for _has_changed
    ok(!$ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives no change');
    
    # create a new lane
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $new_lane = VRTrack::Lane->create($vrtrack, 'new_lane');
    $new_lane->is_withdrawn(0);
    $new_lane->library_id(16);
    $new_lane->update;
    ok($ds->_source_instance->_has_changed, 'vrtrack datasource _has_changed gives change after new lane insertion in test vrtrack db');
    
    # go back to previous state by deleting this lane
    $new_lane->delete;
    $ds->_source_instance->_has_changed;
    ok(!$ds->_source_instance->_has_changed, 'vrtrack datasource _has_change gives no change after inserted lane deleted in test vrtrack db');
    
    # add a file to another lane and check for changes
    my $lane_to_add_file_for = VRTrack::Lane->new_by_name($vrtrack, 'ERR003040');
    my $newfile = $lane_to_add_file_for->add_file('new.fastq');
    is $ds->_source_instance->_has_changed, 0, 'vrtrack datasource _has_changed gives no change after adding a file in test vrtrack db, with method lanes';
    $newfile->delete;
    
    # change some md5 sums in the files
    my $file = VRTrack::File->new_by_hierarchy_name($vrtrack, 'ERR003038.filt.fastq.gz');
    # check for changes
    $file->md5('34c009157187c5d9a7e976563ec1bad8');
    $file->update;
    is $ds->_source_instance->_has_changed, 0, 'datasource _has_changed got no change after md5 change in file table in test vrtrack db, with method lanes';
    $file->md5('cac33e4fc8ff2801978cfd5a223f5064');
    $file->update;
    
    # if we change something that doesn't indicate a real change as far as our
    # options are concerned, it shouldn't come up as _has_changed
    $lane_to_add_file_for->qc_status('pending');
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 0, 'changing qc_status when it was not an option is ignored';
    $lane_to_add_file_for->auto_qc_status('failed');
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 1, 'changing auto_qc_status when it was an option is detected';
    $lane_to_add_file_for->auto_qc_status('no_qc');
    $lane_to_add_file_for->qc_status('no_qc');
    $lane_to_add_file_for->update;
    $ds->_source_instance->_has_changed;
    is $ds->_source_instance->_has_changed, 0, 'reverting lane back gives no change';
    $lane_to_add_file_for->is_withdrawn(1);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 1, 'changing withdrawn is always detected';
    $lane_to_add_file_for->is_withdrawn(0);
    $lane_to_add_file_for->update;
    $ds->_source_instance->_has_changed;
    is $ds->_source_instance->_has_changed, 0, 'prior to changing sample, _has_changed is false';
    my $sample = VRTrack::Sample->new($vrtrack, 1);
    $sample->name('foo');
    $sample->update;
    is $ds->_source_instance->_has_changed, 1, 'changing sample name is always detected';
    $sample->name('NA07056');
    $sample->update;
    my $individual = VRTrack::Individual->new($vrtrack, 1);
    $individual->name('foo');
    $individual->update;
    is $ds->_source_instance->_has_changed, 1, 'changing individual name is always detected';
    $individual->name('NA07056');
    $individual->update;
    $individual->alias('foo');
    $individual->update;
    is $ds->_source_instance->_has_changed, 1, 'changing individual alias is always detected';
    $individual->alias('NA07056');
    $individual->update;
    
    $ds = VRPipe::DataSource->create(
        type    => 'vrtrack',
        method  => 'lanes',
        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
        options => { qc_status => 'investigate' }
    );
    $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 0, 'initially got no results for vrtrack lanes qc_status => investigate';
    $lane_to_add_file_for->qc_status('investigate');
    $lane_to_add_file_for->update;
    $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 1, 'after changing a lane to qc investigate, got 1 dataelement';
    
    $ds = VRPipe::DataSource->create(
        type    => 'vrtrack',
        method  => 'lanes',
        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
        options => {}
    );
    $results = 0;
    foreach my $element (@{ get_elements($ds) }) {
        $results++;
    }
    is $results, 55, 'with no options we get all the active lanes in the db';
    $lane_to_add_file_for->qc_status('pending');
    $lane_to_add_file_for->is_processed('mapped', 1);
    $lane_to_add_file_for->raw_bases(50);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 0, 'we can change lots of things without being detected';
    $lane_to_add_file_for->is_withdrawn(1);
    $lane_to_add_file_for->update;
    is $ds->_source_instance->_has_changed, 1, 'but withdrawn is still detected';
    
    $lane_to_add_file_for->qc_status('no_qc');
    $lane_to_add_file_for->is_processed('mapped', 0);
    $lane_to_add_file_for->raw_bases(1473337566);
    $lane_to_add_file_for->is_withdrawn(0);
    $lane_to_add_file_for->update;
    
    # lane_fastqs tests
    ok $ds = VRPipe::DataSource->create(
        type    => 'vrtrack',
        method  => 'lane_fastqs',
        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
        options => { qc_status => 'no_qc|pending', local_root_dir => dir('t')->absolute->stringify, library_regex => 'g1k-sc-NA19190-YRI-1\|SC\|SRP000542\|NA19190' }
      ),
      'could create a vrtrack datasource';
    
    $individual = VRTrack::Individual->new_by_name($vrtrack, 'NA19190');
    $individual->alias('NA19190_from_alias');
    $individual->update;
    
    @results = ();
    foreach my $element (@{ get_elements($ds) }) {
        push(@results, result_with_inflated_paths($element));
    }
    
    is_deeply $results[0],
      {
        paths => [file(qw(t data NA19190 sequence_read ERR003199.filt.fastq.gz))->absolute, file(qw(t data NA19190 sequence_read ERR003199_1.filt.fastq.gz))->absolute, file(qw(t data NA19190 sequence_read ERR003199_2.filt.fastq.gz))->absolute],
        lane  => 'ERR003199'
      },
      'got correct results for vrtrack lane_fastqs';
    my $vrfile = VRPipe::File->get(path => file(qw(t data NA19190 sequence_read ERR003199.filt.fastq.gz))->absolute);
    my $meta = $vrfile->metadata;
    is_deeply $meta,
      {
        'bases'        => '1696783',
        'withdrawn'    => 0,
        'project'      => 'SRP000542',
        'species'      => 'Homo sapiens',
        'population'   => 'YRI',
        'paired'       => 1,
        'reads'        => '45859',
        'library'      => 'g1k_sc_NA19190_YRI_1',
        'lane_id'      => '101',
        'individual'   => 'NA19190_from_alias',
        'center_name'  => 'SC',
        'sample'       => 'NA19190',
        'platform'     => 'SLX',
        'expected_md5' => 'dfa4364855815d7433c451a87f0520d0',
        'study'        => 'SRP000542',
        'lane'         => 'ERR003199',
        'insert_size'  => 175
      },
      'a VRPipe::File created by vrtrack datasource has the correct metadata';
    
    $newfile = $lane_to_add_file_for->add_file('new.fastq');
    $newfile->type(0);
    $newfile->update;
    is $ds->_source_instance->_has_changed, 1, 'vrtrack datasource _has_changed gives change after adding a fastq file in test vrtrack db, with method lane_fastqs';
    $newfile->delete;
    
    $newfile = $lane_to_add_file_for->add_file('new.bam');
    $newfile->type(5);
    $newfile->update;
    $ds->_source_instance->_has_changed;
    is $ds->_source_instance->_has_changed, 0, 'vrtrack datasource _has_changed gives no change after adding a bam file in test vrtrack db, with method lane_fastqs';
    $newfile->delete;
    
    # change some md5 sums in the files
    $file = VRTrack::File->new_by_hierarchy_name($vrtrack, 'ERR003038.filt.fastq.gz');
    # check for changes
    $file->md5('34c009157187c5d9a7e976563ec1bad8');
    $file->update;
    is $ds->_source_instance->_has_changed, 1, 'datasource _has_changed got change after md5 change in file table in test vrtrack db, with method lane_fastqs';
    $file->md5('cac33e4fc8ff2801978cfd5a223f5064');
    $file->update;
    $ds->_source_instance->_has_changed;
    is $ds->_source_instance->_has_changed, 0, 'reverting file md5 back gives no change';
    
    # if we change the insert_size in vrtrack, this should cause the datasource
    # to reset the elements. Critically, the vrfile metadata should get updated,
    # or else we'd be stuck in an infinite loop of always detecting the
    # insert_size change and resetting
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my $lib = VRTrack::Library->new_by_name($vrtrack, 'g1k-sc-NA19190-YRI-1|SC|SRP000542|NA19190');
    $lib->fragment_size_from(200);
    $lib->fragment_size_to(200);
    $lib->update;
    # (we also add a lane to force datasource to notice a change)
    VRTrack::Lane->create($vrtrack, 'new_lane2');
    get_elements($ds);
    $vrfile->reselect_values_from_db;
    is $vrfile->metadata->{insert_size}, '200', 'changing insert_size in vrtrack changes insert_size metadata on vrpipe files';
    
    # test that changing the md5 in the vrtrack db results in the corresponding
    # dataelement being started from scratch, and the expected_md5 on the file
    # getting updated
    my ($test_de) = $vrfile->input_to;
    my $test_des = VRPipe::DataElementState->create(dataelement => $test_de->id, pipelinesetup => 1, completed_steps => 1);
    is $test_des->completed_steps, 1, 'got and pretended a dataelementstate completed a step';
    $file = VRTrack::File->new_by_hierarchy_name($vrtrack, 'ERR003199.filt.fastq.gz');
    $file->md5('foo');
    $file->update;
    get_elements($ds);
    $test_des->reselect_values_from_db;
    is $test_des->completed_steps, 0, 'changing md5 on a dataelement input file started it from scratch';
    $vrfile->reselect_values_from_db;
    is $vrfile->meta_value('expected_md5'), 'foo', 'and it changed the metadata on the file';
    
    # test getting improved bams that passed qc
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    my %expected_qc_passed_improved_bams;
    my %passed_hnames = map { $_ => 1 } @{ $expectations{qc_status_passed} };
    foreach my $hname (@{ $expectations{gt_status_confirmed} }) {
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $hname);
        my $improved_bam = VRPipe::File->create(path => file($hname . '.improved.bam')->absolute);
        my $vrfile_name  = 'VRPipe::File::' . $improved_bam->id;
        my $md5          = 'an_md5_' . $improved_bam->id;
        $vrfile = $vrlane->add_file($vrfile_name);
        $vrfile->type(5);
        $vrfile->md5($md5);
        $vrfile->update;
        $improved_bam->md5($md5);
        $improved_bam->update;
        
        if (exists $passed_hnames{$hname}) {
            $expected_qc_passed_improved_bams{ $improved_bam->path->stringify } = 1;
        }
        
        # to give us more than 1 group to test, change the library of one of the
        # lanes
        if ($hname eq 'ERR008838') {
            $vrlane->library_id(1);
            $vrlane->update;
        }
    }
    my %expected_groups = (
        'SRP000546|SRP000546|NA18633|HUMsgR3AIDCAASE'      => 1,
        'SRP000547|SRP000547|NA07056|g1k_sc_NA07056_CEU_1' => 1
    );
    $ds = VRPipe::DataSource->create(
        type    => 'vrtrack',
        method  => 'lane_improved_bams',
        source  => $ENV{VRPIPE_VRTRACK_TESTDB},
        options => {
            gt_status         => 'confirmed',
            qc_status         => 'passed',
            group_by_metadata => 'project|study|sample|library'
        }
    );
    my %actual_qc_passed_improved_bams;
    my %actual_groups;
    foreach my $element (@{ get_elements($ds) }) {
        my $result = result_with_inflated_paths($element);
        $actual_groups{ $result->{group} || 'no_group' } = 1;
        foreach my $path (@{ $result->{paths} }) {
            $actual_qc_passed_improved_bams{$path} = 1;
        }
    }
    is_deeply \%actual_qc_passed_improved_bams, \%expected_qc_passed_improved_bams, 'an improved bam qc passed datasource gave the expected files';
    is_deeply \%actual_groups, \%expected_groups, 'a group_by_metadata datasource gave the expected groups';
}

exit;
