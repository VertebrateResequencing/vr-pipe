#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => 'VRPIPE_TEST_PIPELINES',
		    required_exe => [qw(glf samtools bin2hapmap)]);
    use TestPipelines;
}

my $checking_output_dir = get_output_dir('bam_genotype_checking');

ok my $gc_pipeline = VRPipe::Pipeline->create(name => 'bam_genotype_checking'), 'able to get a pre-written pipeline';

my @s_names;
foreach my $stepmember ($gc_pipeline->steps) {
    push(@s_names, $stepmember->step->name);
}
my @expected_step_names = qw(bin2hapmap_sites
                             mpileup_bcf_hapmap
			     glf_check_genotype
			     gtypex_genotype_analysis);

is_deeply \@s_names, \@expected_step_names, 'the pipeline has the correct steps';

my $snp_bin_source = file(qw(t data qc_chr20_snps.bin));
my $snp_dir = dir($checking_output_dir, 'snp');
$gc_pipeline->make_path($snp_dir);
my $snp_bin = file($snp_dir, 'qc_chr20_snps.bin')->stringify;
copy($snp_bin_source, $snp_bin);

my $ref_fa_source = file(qw(t data human_g1k_v37.chr20.fa));
my $ref_dir = dir($checking_output_dir, 'ref');
$gc_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'human_g1k_v37.chr20.fa')->stringify;
copy($ref_fa_source, $ref_fa);

my $ds = VRPipe::DataSource->create(type => 'fofn',
			         method => 'all',
				 source => file(qw(t data hs_chr20.qc.bam.fofn))->absolute);

VRPipe::PipelineSetup->create(name => 'genotype_checking',
			   datasource => $ds,
			   output_root => $checking_output_dir,
			   pipeline => $gc_pipeline,
			   options => {reference_fasta => $ref_fa,
			               hapmap2bin_sample_genotypes_file => $snp_bin,
				       expected_sample_metadata_key => 'sample'});

my (@output_files,@final_files);
my $element_id=0;
foreach my $in ('hs_chr20.a', 'hs_chr20.b', 'hs_chr20.c', 'hs_chr20.d') {
    # mpileup_bcf_hapmap step will add sample metadata by figuring it out from
    # bam; we'll add individual metadata to test glf_check_gentotype's
    # expected_sample_metadata_key option, and to test we get confirmed results
    # sometimes
    VRPipe::File->create(path => file('t', 'data', $in.'.bam')->absolute)->add_metadata({individual => 'NA20586'});
    
    $element_id++;
    my @output_subdirs = output_subdirs($element_id);
    push(@output_files, file(@output_subdirs, '2_mpileup_bcf_hapmap', "${in}.bam.bcf"));
    push(@final_files, file(@output_subdirs, '3_glf_check_genotype', "${in}.bam.bcf.gtypex"));
}

ok handle_pipeline(@output_files, @final_files), 'pipeline ran and created all expected output files';

my @gtype_results;
foreach (qw(a b c d)) {
    my $bam_meta = VRPipe::File->create(path => file('t', 'data', "hs_chr20.$_.bam")->absolute)->metadata;
    push(@gtype_results, $bam_meta->{gtype_analysis});
}
is_deeply \@gtype_results, ['status=unconfirmed expected=NA20526 found=NA20586 ratio=1.016',
			    'status=unconfirmed expected=NA20527 found=NA20521 ratio=1.000',
			    'status=unconfirmed expected=NA20526 found=NA20521 ratio=1.000',
			    'status=unconfirmed expected=NA20526 found=NA20521 ratio=1.000'], 'gtype_analysis results were stored correctly as metadata on the input bams';

VRPipe::PipelineSetup->create(name => 'genotype_checking',
			   datasource => $ds,
			   output_root => $checking_output_dir,
			   pipeline => $gc_pipeline,
			   options => {reference_fasta => $ref_fa,
			               hapmap2bin_sample_genotypes_file => $snp_bin,
				       gtype_confidence => 1.01}); # expected_sample_metadata_key defaults to individual

ok handle_pipeline(), 'pipeline ran again';

@gtype_results = ();
foreach (qw(a b c d)) {
    my $bam_meta = VRPipe::File->create(path => file('t', 'data', "hs_chr20.$_.bam")->absolute)->metadata;
    push(@gtype_results, $bam_meta->{gtype_analysis});
}
is_deeply \@gtype_results, ['status=confirmed expected=NA20586 found=NA20586 ratio=1.016',
			    'status=unconfirmed expected=NA20586 found=NA20521 ratio=1.000',
			    'status=unconfirmed expected=NA20586 found=NA20521 ratio=1.000',
			    'status=unconfirmed expected=NA20586 found=NA20521 ratio=1.000'], 'gtype_analysis results were updated correctly when we reran using individual as the expected';

done_testing;

finish;