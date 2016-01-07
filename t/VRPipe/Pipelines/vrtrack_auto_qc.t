#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;
use File::Copy;
use VRPipe::Parser;

BEGIN {
    use Test::Most tests => 16;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES VRPIPE_VRTRACK_TESTDB)],
        required_exe => [qw(samtools bcftools bamcheck)]                  # latest versions required
    );
    use TestPipelines;
    
    use_ok('VRTrack::Factory');
    use_ok('VRPipe::Steps::bamcheck');
    use_ok('VRPipe::Steps::vrtrack_update_mapstats');
}

# setup a little VRTrack db to hold info on our test bams
my @lane_names = ('7046_8#8', '7936_7#52');
my %lib_names = ('7046_8#8' => 'KUU25220302_3887679', '7936_7#52' => '5002874');
my $output_dir = get_output_dir('auto_qc');
my $vrtrack;
prepare_vrtrack();

my $auto_qc_tab_file = file(qw(t data autoqc.tab))->absolute->stringify;

# autoqc that writes results to vrtrack
VRPipe::PipelineSetup->create(
    name       => 'auto_qc',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        source  => $auto_qc_tab_file,
        options => { delimiter => "\t" }
    ),
    output_root => $output_dir,
    pipeline    => VRPipe::Pipeline->create(name => 'vrtrack_auto_qc'),
    options     => { vrtrack_db => $ENV{VRPIPE_VRTRACK_TESTDB}, auto_qc_min_ins_to_del_ratio => 0.5 }
);

ok handle_pipeline(), 'vrtrack_auto_qc pipeline ran';

# check autoqc tests on vrtrack database
my %actual_auto_qc_data;
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
die "Can't connect to tracking database\n" unless $vrtrack;

foreach my $lane_name (@lane_names) {
    my $sql = qq[select test,CASE result when 1 then "PASSED" when 0 then "FAILED" END as result,reason,m.mapstats_id from autoqc as a join latest_mapstats as m on m.mapstats_id = a.mapstats_id join latest_lane as l on l.lane_id = m.lane_id where l.name='$lane_name' order by a.autoqc_id;];
    
    my $sth = $vrtrack->{_dbh}->prepare($sql);
    $sth->execute();
    while (my $r = $sth->fetchrow_hashref) {
        push @{ $actual_auto_qc_data{$lane_name} }, "$r->{test}:\t$r->{result}\t # $r->{reason}\n";
    }
    $sth->finish;
}

my %expected_auto_qc_data;
foreach my $lane_name (@lane_names) {
    my $file_name = $lane_name;
    $file_name =~ s/#/_/;
    my $aqcfile = file('t', 'data', $file_name . '.auto_qc.txt')->absolute;
    open FILE, "<", $aqcfile or die $aqcfile, $!;
    while (<FILE>) {
        push @{ $expected_auto_qc_data{$lane_name} }, $_ unless /^Verdict/;
    }
    close FILE;
}
is_deeply \%actual_auto_qc_data, \%expected_auto_qc_data, 'auto qc pipeline generated the expected autoqc on vrtrack showing why the lanes passed';

my $passed_auto_qc_lanes = 0;
my $failed_auto_qc_libs  = 0;
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
foreach my $lane_name (@lane_names) {
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
    $passed_auto_qc_lanes++ if $lane->auto_qc_status eq 'passed';
    my $lib = VRTrack::Library->new($vrtrack, $lane->library_id);
    $failed_auto_qc_libs++ if $lib->auto_qc_status eq 'failed';
}
is_deeply [$passed_auto_qc_lanes, $failed_auto_qc_libs], [0, 2], 'auto qc pipeline set lanes and libraries in VRTrack to passed/failed as appropriate';

# also test the combined genotype check + auto qc pipeline
$output_dir = get_output_dir('genotype_and_auto_qc');

my $aqg_pipeline = VRPipe::Pipeline->create(name => 'vrtrack_auto_qc_with_genotype_checking');
my @s_names;
foreach my $stepmember ($aqg_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(bam_index bcftools_generate_sites_file mpileup_vcf vcf_index bcftools_gtcheck bcftools_genotype_analysis vrtrack_auto_qc);
is_deeply \@s_names, \@expected_step_names, 'the vrtrack_auto_qc_with_genotype_checking pipeline has the correct steps';

my $geno_source = file(qw(t data hs_chr1.genotypes.bcf));
my $geno_dir = dir($output_dir, 'geno');
$aqg_pipeline->make_path($geno_dir);
my $geno = file($geno_dir, 'hs_chr1.genotypes.bcf')->stringify;
copy($geno_source,       $geno);
copy("$geno_source.csi", "$geno.csi");

my $tab_file_with_geno   = file($geno_dir, 'delimited_datasource.txt');
my $tab_file_with_geno_2 = file($geno_dir, 'delimited_datasource_2.txt');
my $tfwg_file  = VRPipe::File->create(path => $tab_file_with_geno);
my $tfwg2_file = VRPipe::File->create(path => $tab_file_with_geno_2);
my $tfwg_ofh   = $tfwg_file->openw;
my $tfwg2_ofh  = $tfwg2_file->openw;
my $aqtf_file  = VRPipe::File->create(path => $auto_qc_tab_file);
my $aqtf_ifh   = $aqtf_file->openr;
my @public_names = qw(foo bar);

while (<$aqtf_ifh>) {
    chomp;
    print $tfwg_ofh $_, "\t$geno\n";
    $_ =~ s/bamcheck/bai/;
    print $tfwg2_ofh $_, "\t$geno\n";
    my ($bam) = split;
    my $file = VRPipe::File->get(path => file($bam)->absolute);
    $file->add_metadata({ public_name => pop(@public_names) });
}
$aqtf_file->close;
$tfwg_file->close;
$tfwg2_file->close;

my $ref_fa_source = file(qw(t data human.chr1_chunk.fa));
my $ref_dir = dir($output_dir, 'ref');
$aqg_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human.chr1_chunk.fa')->stringify;
copy($ref_fa_source, $ref_fa);

prepare_vrtrack();

VRPipe::PipelineSetup->create(
    name       => 'auto_qc with genotype checking',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        source  => $tab_file_with_geno,
        options => { delimiter => "\t" }
    ),
    output_root => $output_dir,
    pipeline    => $aqg_pipeline,
    options     => {
        vrtrack_db                        => $ENV{VRPIPE_VRTRACK_TESTDB},
        auto_qc_min_ins_to_del_ratio      => 0.5,
        reference_fasta                   => $ref_fa,
        vcf_sample_from_metadata          => 'sample:public_name+sample',
        expected_sample_from_metadata_key => 'public_name+sample',
        cleanup                           => 0
    }
);

my (@output_files, @final_files);
my $element_id = 2;
foreach my $in ('autoqc_normal', 'autoqc_short') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 2);
    push(@output_files, file(@output_subdirs, '3_mpileup_vcf',      "mpileup.vcf.gz"));
    push(@output_files, file(@output_subdirs, '3_mpileup_vcf',      "mpileup.vcf.gz"));
    push(@final_files,  file(@output_subdirs, '5_bcftools_gtcheck', "mpileup.vcf.gz.gtypex"));
}

ok handle_pipeline(@output_files, @final_files), 'vrtrack_auto_qc_with_genotype_checking pipeline ran and created all expected output files';

my @found_gtype_analysis;
foreach my $in ('autoqc_normal', 'autoqc_short') {
    my $meta = VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->metadata;
    push(@found_gtype_analysis, $meta->{gtype_analysis});

}
my @expected_gtype_analysis = (
    'status=confirmed expected=bar_KUU25220302 found=bar_KUU25220302 ratio=1.00 concordance=1.000',
    'status=confirmed expected=foo_SISu5277216 found=foo_SISu5277216 ratio=1.00 concordance=1.000',
);

is_deeply \@found_gtype_analysis, \@expected_gtype_analysis, 'pipeline generated the expected genotype analysis metadata';

# check autoqc tests on vrtrack database
%actual_auto_qc_data = ();
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
die "Can't connect to tracking database\n" unless $vrtrack;

foreach my $lane_name (@lane_names) {
    my $sql = qq[select test,CASE result when 1 then "PASSED" when 0 then "FAILED" END as result,reason,m.mapstats_id from autoqc as a join latest_mapstats as m on m.mapstats_id = a.mapstats_id join latest_lane as l on l.lane_id = m.lane_id where l.name='$lane_name' order by a.autoqc_id;];
    
    my $sth = $vrtrack->{_dbh}->prepare($sql);
    $sth->execute();
    while (my $r = $sth->fetchrow_hashref) {
        push @{ $actual_auto_qc_data{$lane_name} }, "$r->{test}:\t$r->{result}\t # $r->{reason}\n";
    }
    $sth->finish;
    
    unshift(@{ $expected_auto_qc_data{$lane_name} }, "Genotype check:\tPASSED\t # The status is 'confirmed'.\n");
}

is_deeply \%actual_auto_qc_data, \%expected_auto_qc_data, 'auto qc pipeline generated the expected autoqc on vrtrack showing why the lanes passed';

# also test the combined genotype check + auto qc pipeline based on summary statistics stored in the graph db
$output_dir = get_output_dir('genotype_and_graph_auto_qc');

$aqg_pipeline = VRPipe::Pipeline->create(name => 'vrtrack_auto_qc_from_graph_db_with_genotype_checking');
@s_names = ();
foreach my $stepmember ($aqg_pipeline->step_members) {
    push(@s_names, $stepmember->step->name);
}
@expected_step_names = qw(bcftools_generate_sites_file mpileup_vcf bcftools_gtcheck bcftools_genotype_analysis graph_auto_qc update_vrtrack_from_graph_db);
is_deeply \@s_names, \@expected_step_names, 'the vrtrack_auto_qc_from_graph_db_with_genotype_checking pipeline has the correct steps';

# parse the bamcheck files and add the summary statistics to Bam_Stats nodes of input bam files along with some extra information that graph_auto_qc step requires and are normally generated in advance by npg_cram_stats_parser step.
ok my $schema = VRPipe::Schema->create('VRTrack'), 'able to create a VRTrack schema instance';
my $vrpipe = VRPipe::Schema->create('VRPipe');
foreach my $in ('autoqc_normal', 'autoqc_short') {
    my $bam_file = $vrpipe->add('File', { path => file('t', 'data', $in . '.bam')->absolute });
    my $qc_file  = $vrpipe->add('File', { path => file('t', 'data', $in . '.bam.bamcheck')->absolute });
    my $bc = VRPipe::Parser->create('bamcheck', { file => file('t', 'data', $in . '.bam.bamcheck')->absolute });
    my %stats = ();
    foreach my $key (keys %$bc) {
        my $nk = $key;
        $nk =~ s/\_/ /g;
        $nk =~ s/cigar/\(cigar\)/;
        my $val = $bc->$key;
        if (ref($val) eq 'ARRAY') {
            my @vals;
            foreach my $v (@$val) {
                push(@vals, "$v");
            }
            $val = \@vals;
        }
        elsif (ref($val) eq 'HASH') {
            next;
        }
        $stats{$nk} = $val;
    }
    $bam_file->relate_to($qc_file, 'qc_file');
    $schema->add(
        'Bam_Stats',
        {
            mode                                    => 'rmdup',
            options                                 => '-opt',
            date                                    => time(),
            'median of indel cycle counts'          => '3,1,3,4',                                   #fake values
            'max of indel cycle counts'             => '11,10,9,14',                                #fake values
            'dupl mapped bases'                     => ($in eq 'autoqc_normal') ? '570' : '5411',
            'total mapped bases'                    => ($in eq 'autoqc_normal') ? '81800' : '84800',
            'reads mapped after rmdup'              => $bc->reads_mapped - $bc->reads_duplicated,
            'inserts within 25% max peak'           => ($in eq 'autoqc_normal') ? '76.04' : '37.97',
            'peak window containing 80% of inserts' => ($in eq 'autoqc_normal') ? '27.73' : '69.44',
            'number of insertions'                  => ($in eq 'autoqc_normal') ? '11' : '28',
            'number of deletions'                   => ($in eq 'autoqc_normal') ? '17' : '43',
            'lane_is_paired_read'                   => 1,
            %stats,
        },
        incoming => { type => 'summary_stats', node => $qc_file }
    );
}

prepare_vrtrack();

VRPipe::PipelineSetup->create(
    name       => 'graph db auto_qc and update vrtrack with genotype checking',
    datasource => VRPipe::DataSource->create(
        type    => 'delimited',
        method  => 'all_columns',
        source  => $tab_file_with_geno_2,
        options => { delimiter => "\t" }
    ),
    output_root => $output_dir,
    pipeline    => $aqg_pipeline,
    options     => {
        vrtrack_db                        => $ENV{VRPIPE_VRTRACK_TESTDB},
        auto_qc_min_ins_to_del_ratio      => 0.5,
        reference_fasta                   => $ref_fa,
        vcf_sample_from_metadata          => 'sample:public_name+sample',
        expected_sample_from_metadata_key => 'public_name+sample',
        cleanup                           => 0,
    }
);

ok handle_pipeline(), 'vrtrack_auto_qc_from_graph_db_with_genotype_checking pipeline ran';

@output_files = ();
@final_files  = ();
$element_id   = 4;
foreach my $in ('autoqc_normal', 'autoqc_short') {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id, 3);
    push(@output_files, file(@output_subdirs, '2_mpileup_vcf',      "mpileup.vcf.gz"));
    push(@final_files,  file(@output_subdirs, '3_bcftools_gtcheck', "mpileup.vcf.gz.gtypex"));
}

ok handle_pipeline(@output_files, @final_files), 'vrtrack_auto_qc_from_graph_db_with_genotype_checking pipeline ran and created all expected output files';

@found_gtype_analysis = ();
foreach my $in ('autoqc_normal', 'autoqc_short') {
    my $meta = VRPipe::File->create(path => file('t', 'data', $in . '.bam')->absolute)->metadata;
    push(@found_gtype_analysis, $meta->{gtype_analysis});

}

is_deeply \@found_gtype_analysis, \@expected_gtype_analysis, 'pipeline generated the expected genotype analysis metadata';

# check autoqc tests on vrtrack database
%actual_auto_qc_data = ();
$vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'r');
die "Can't connect to tracking database\n" unless $vrtrack;

foreach my $lane_name (@lane_names) {
    my $sql = qq[select test,CASE result when 1 then "PASSED" when 0 then "FAILED" END as result,reason,m.mapstats_id from autoqc as a join latest_mapstats as m on m.mapstats_id = a.mapstats_id join latest_lane as l on l.lane_id = m.lane_id where l.name='$lane_name' order by a.autoqc_id;];
    
    my $sth = $vrtrack->{_dbh}->prepare($sql);
    $sth->execute();
    while (my $r = $sth->fetchrow_hashref) {
        push @{ $actual_auto_qc_data{$lane_name} }, "$r->{test}:\t$r->{result}\t # $r->{reason}\n";
    }
    $sth->finish;
    
    @{ $actual_auto_qc_data{$lane_name} }   = sort { $a cmp $b } @{ $actual_auto_qc_data{$lane_name} };
    @{ $expected_auto_qc_data{$lane_name} } = sort { $a cmp $b } @{ $expected_auto_qc_data{$lane_name} };
    splice(@{ $expected_auto_qc_data{$lane_name} }, 9, 0, "NPG QC status check:\tFAILED\t # The lane failed the NPG QC check, so we auto-fail as well since this data will not be auto-submitted to EGA/ENA.\n");
}

is_deeply \%actual_auto_qc_data, \%expected_auto_qc_data, 'auto qc pipeline generated the expected autoqc on vrtrack showing why the lanes passed';

finish;

sub prepare_vrtrack {
    my %cd = VRTrack::Factory->connection_details('rw');
    open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
    print $mysqlfh "drop database if exists $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "create database $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    print $mysqlfh "use $ENV{VRPIPE_VRTRACK_TESTDB};\n";
    my @sql = VRTrack::VRTrack->schema();
    foreach my $sql (@sql) {
        print $mysqlfh $sql;
    }
    close($mysqlfh);
    
    $vrtrack = VRTrack::Factory->instantiate(database => $ENV{VRPIPE_VRTRACK_TESTDB}, mode => 'rw');
    
    foreach my $lane_name (@lane_names) {
        my $lane = VRTrack::Lane->create($vrtrack, $lane_name);
        $lane->hierarchy_name($lane_name);
        my $lib = VRTrack::Library->create($vrtrack, $lib_names{$lane_name});
        $lane->library_id($lib->id);
        $lane->update;
    }
    
    # vrtrack_auto_qc step requires input files have lane metadata; it also requires
    # that it have a populated mapstats table, which in turn means we need the
    # bamcheck results stored as metadata on the bams
    my $lni = 0;
    foreach my $file_prefix (qw(autoqc_normal autoqc_short)) {
        foreach my $suffix ('.bam', '.bam.bamcheck') {
            my $file = VRPipe::File->create(path => file('t', 'data', $file_prefix . $suffix)->absolute);
            $file->add_metadata({ lane => $lane_names[$lni] });
        }
        
        # bamcheck metadata
        my $bam = VRPipe::File->create(path => file('t', 'data', $file_prefix . '.bam')->absolute)->path;
        my $bamcheck = VRPipe::File->create(path => file($output_dir, $file_prefix . '.bam.bamcheck'))->path;
        VRPipe::Steps::bamcheck->stats_from_bamcheck("bamcheck $bam > $bamcheck");
        
        # mapstats
        VRPipe::Steps::vrtrack_update_mapstats->update_mapstats(db => $ENV{VRPIPE_VRTRACK_TESTDB}, bam => $bam, lane => $lane_names[$lni], plot_dir => $output_dir, plots => []);
        
        $lni++;
    }
}
