#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use POSIX qw(getgroups);

BEGIN {
    use Test::Most tests => 13;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_to_fastq');
    use_ok('VRPipe::Parser');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_to_fastq', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Converts bam files to fastq files'], 'bam_to_fastq step created and has correct description';

# pipeline requires certain metadata on the input bams, so we just manually set that now
VRPipe::File->create(
    path     => file(qw(t data 2822_6.pe.bam))->absolute,
    metadata => {
        lane          => '2822_6',
        reads         => 400,
        bases         => 23000,
        forward_reads => 200,
        reverse_reads => 200,
        paired        => 1
    }
);
VRPipe::File->create(
    path     => file(qw(t data 2822_6.improved.pe.bam))->absolute,
    metadata => {
        lane          => '2822_6',
        reads         => 400,
        bases         => 23000,
        forward_reads => 200,
        reverse_reads => 200,
        paired        => 1
    }
);
VRPipe::File->create(
    path     => file(qw(t data NA19381.bam))->absolute,
    metadata => {
        sample        => 'NA19381',
        reads         => 80,
        bases         => 8080,
        forward_reads => 40,
        reverse_reads => 40,
        paired        => 1
    }
);
VRPipe::File->create(
    path     => file(qw(t data 2822_6.se.pe.bam))->absolute,
    metadata => {
        lane          => '2822_6',
        reads         => 400,
        bases         => 23000,
        forward_reads => 200,
        reverse_reads => 200,
        paired        => 0
    }
);

# since this step produces a temp file, we'll also test that setting unix_group
# option on the pipelinesetup doesn't break everything
my @groups = map { scalar getgrgid($_) } getgroups();

my $setup = VRPipe::PipelineSetup->create(
    name        => 'btq_setup',
    datasource  => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data datasource.btf_fofn))->absolute),
    output_root => $output_dir,
    unix_group  => $groups[0],
    pipeline    => $pipeline,
    options     => {}
);

my @fastqs;
my %bases = (1 => '2822_6.pe_', 2 => '2822_6.improved.pe_', 3 => 'NA19381_', 4 => '2822_6.se.pe_');
for my $j (1, 2, 3, 4) {
    foreach my $i (1, 2) {
        my @output_subdirs = output_subdirs($j);
        push(@fastqs, file(@output_subdirs, '1_bam_to_fastq', "$bases{$j}$i.fastq"));
    }
}
push(@fastqs, file(output_subdirs(4), '1_bam_to_fastq', "2822_6.se.pe_M.fastq"));
ok handle_pipeline(@fastqs), 'bam_to_fastq pipeline ran ok, generating the expected fastqs';

is_deeply read_fastq(file(qw(t data 2822_6_1.fastq))->absolute), read_fastq($fastqs[0]), 'forward fastq had same data as original forward fastq the input bam was mapped from';
is_deeply read_fastq(file(qw(t data 2822_6_2.fastq))->absolute), read_fastq($fastqs[1]), 'reverse fastq had same data as original reverse fastq the input bam was mapped from';
is_deeply read_fastq(file(qw(t data 2822_6_1.fastq))->absolute), read_fastq($fastqs[2]), 'forward fastq from improved bam had same data as original forward fastq the input bam was mapped from';
is_deeply read_fastq(file(qw(t data 2822_6_2.fastq))->absolute), read_fastq($fastqs[3]), 'reverse fastq from improved bam had same data as original reverse fastq the input bam was mapped from';
is_deeply [scalar(keys %{ read_fastq($fastqs[4]) }), scalar(keys %{ read_fastq($fastqs[5]) }), VRPipe::File->get(path => $fastqs[4])->meta_value('sample')], [40, 40, 'NA19381'], 'merged bam produced correctly sized fastqs with sample metadata attached';

my $fmeta = VRPipe::File->create(path => $fastqs[0])->metadata;
is_deeply [$fmeta->{reads}, $fmeta->{bases}, $fmeta->{avg_read_length}], [200, 12200, '61.00'], 'basic metadata is correct for the forward fastq';

my %reads = (1 => 199, 2 => 199, 'M' => 2);
for my $i (1, 2, 'M') {
    my $path = file(output_subdirs(4), '1_bam_to_fastq', "2822_6.se.pe_$i.fastq")->stringify;
    my $f = VRPipe::File->get(path => "$path");
    is_deeply [$reads{$i}], [$f->metadata->{reads}], "$i fastq generated with correct reads metadata for a bam with both SE and PE reads";
}

finish;

sub read_fastq {
    my $path = shift;
    my $pars = VRPipe::Parser->create('fastq', { file => VRPipe::File->create(path => $path) });
    
    my %data;
    my $pr = $pars->parsed_record();
    while ($pars->next_record()) {
        $data{ $pr->[0] } = [$pr->[1], $pr->[2]];
    }
    return \%data;
}
