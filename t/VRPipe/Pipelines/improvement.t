#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES GATK PICARD)],
                    required_exe => [qw(samtools)]);
    use TestPipelines;
}

my $improvement_output_dir = get_output_dir('bam_improvement');

ok my $improvement_pipeline = VRPipe::Pipeline->create(name => 'bam_improvement'), 'able to get the bam_improvement pipeline';

my @s_names;
foreach my $stepmember ($improvement_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(sequence_dictionary
                             bam_metadata
                             bam_index
                             gatk_target_interval_creator
                             bam_realignment_around_known_indels
                             bam_index
                             bam_count_covariates
                             bam_recalibrate_quality_scores
                             bam_calculate_bq
                             bam_reheader);
is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $ref_fa_source = file(qw(t data S_suis_P17.fa));
my $ref_dir = dir($improvement_output_dir, 'ref');
$improvement_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'S_suis_P17.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $res_dir = dir($improvement_output_dir, 'resources');
$improvement_pipeline->make_path($res_dir);

my $known_indels_source = file(qw(t data known_indels.vcf.gz));
my $known_indels = file($res_dir, 'known_indels.vcf.gz')->stringify;
copy($known_indels_source, $known_indels);
copy($known_indels_source.'.tbi', $known_indels.'.tbi');

my $known_sites_source = file(qw(t data known_sites.vcf.gz));
my $known_sites = file($res_dir, 'known_sites.vcf.gz')->stringify;
copy($known_sites_source, $known_sites);
copy($known_sites_source.'.tbi', $known_sites.'.tbi');

# copy input bams to the output dir, since we will create .bai files and don't
# want them in the t/data directory
my $orig_fofn_file = VRPipe::File->create(path => file(qw(t data improvement_datasource.fofn))->absolute);
my $fofn_file = VRPipe::File->create(path => file($improvement_output_dir, 'improvement_datasource.fofn'));
my $ifh = $orig_fofn_file->openr;
my $ofh = $fofn_file->openw;
while (<$ifh>) {
    chomp;
    my $source = file($_);
    my $dest = file($improvement_output_dir, $source->basename);
    copy($source, $dest);
    print $ofh $dest, "\n";
}
$orig_fofn_file->close;
$fofn_file->close;

ok my $ds = VRPipe::DataSource->create(type => 'fofn',
                                 method => 'all',
                                 source => $fofn_file->path->stringify,
                                 options => {}), 'could create a fofn datasource';

my @results = ();
foreach my $element (@{get_elements($ds)}) {
    push(@results, $element->result);
}
is_deeply \@results, [{paths => [file($improvement_output_dir, '2822_7.pe.bam')->absolute]}, 
                      {paths => [file($improvement_output_dir, '2822_6.pe.bam')->absolute]}, 
                      {paths => [file($improvement_output_dir, '2822_6.se.bam')->absolute]}, 
                      {paths => [file($improvement_output_dir, '2823_4.pe.bam')->absolute]}, 
                      {paths => [file($improvement_output_dir, '8324_8.pe.bam')->absolute]}], 'got correct results for fofn all';

my $improvement_pipelinesetup = VRPipe::PipelineSetup->create(name => 's_suis improvement',
                                                       datasource => $ds,
                                                       output_root => $improvement_output_dir,
                                                       pipeline => $improvement_pipeline,
                                                       options => {reference_fasta => $ref_fa,
                                                                   reference_assembly_name => 'SSuis1',
                                                                   reference_public_url => 'ftp://s.suis.com/ref.fa',
                                                                   reference_species => 'S.Suis',
                                                                   known_indels_for_realignment => "-known $known_indels",
                                                                   known_sites_for_recalibration => "-knownSites $known_sites",
                                                                   gatk_count_covariates_options => '-l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate',
                                                                   gatk_path => $ENV{GATK},
                                                                   picard_path => $ENV{PICARD},
                                                                   cleanup => 0,
                                                                   sequence_dictionary_memory => 150,
                                                                   sequence_dictionary_time => 1});


my @output_files = (file($improvement_output_dir, 'ref', 'S_suis_P17.fa.dict'), file($improvement_output_dir, 'resources', 'known_indels.intervals'));

my @files = ('2822_7.pe.bam', '2822_6.pe.bam', '2822_6.se.bam', '2823_4.pe.bam', '8324_8.pe.bam');
my $element_id = 0;
foreach my $file (@files) {
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    
    push @output_files, file($improvement_output_dir, $file);
    push @output_files, file($improvement_output_dir, "$file.bai");
    $file =~ s/bam$/realign.bam/;
    push @output_files, file(@output_subdirs, '5_bam_realignment_around_known_indels', $file);
    $file =~ s/bam$/recal_data.csv/;
    push @output_files, file(@output_subdirs, '7_bam_count_covariates', $file);
    $file =~ s/recal\_data\.csv$/recal.bam/;
    push @output_files, file(@output_subdirs, '8_bam_recalibrate_quality_scores', $file);
    $file =~ s/bam$/calmd.bam/;
    push @output_files, file(@output_subdirs, '9_bam_calculate_bq', $file);
    push @output_files, file(@output_subdirs, '10_bam_reheader', $file);
}
ok handle_pipeline(@output_files), 'pipeline ran and created all expected output files';

verify_bam_header($output_files[-1], VRPipe::File->create(path => file(qw(t data 8324_8.pe.realign.recal.calmd.bam.reheader))->absolute));

finish;

sub verify_bam_header {
    my ($bam_path, $reheader_file) = @_;
    my %expected_lines;
    parse_header_lines($reheader_file->openr, \%expected_lines);
    
    open(my $sfh, "samtools view -H $bam_path |") || die "Could not view $bam_path with samtools\n";
    my %actual_lines;
    parse_header_lines($sfh, \%actual_lines);
    
    is_deeply \%actual_lines, \%expected_lines, 'header of bam matched expectation';
}

sub parse_header_lines {
    my ($fh, $hash) = @_;
    
    while (<$fh>) {
        chomp;
        next unless /\S/;
        my ($type, $content) = /^\@(\w\w)\t(.+)/;
        
        if ($type eq 'PG') {
            $content =~ s/VN:[^\t]+/VN:version/;
        }
        
        push(@{$hash->{$type}}, $content);
    }
    close($fh);
}