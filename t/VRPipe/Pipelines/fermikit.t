#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use Path::Class;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest (
        required_env => [qw(VRPIPE_TEST_PIPELINES)],
        required_exe => [qw(samtools bfc ropebwt2 fermi2 htsbox bwa vcf-sort)]
    );
    use TestPipelines;
}

ok my $assembly_pipeline = VRPipe::Pipeline->create(name => 'fermikit_unitig_assembly'), 'able to get the fermikit_unitig_assembly pipeline';
ok my $calling_pipeline  = VRPipe::Pipeline->create(name => 'fermikit_calling'),         'able to get the fermikit_calling pipeline';

my $output_dir = get_output_dir('fermikit_test');

my $ref_fa_source = file(qw(t data pombe_ref.fa));
my $ref_dir = dir($output_dir, 'ref');
$assembly_pipeline->make_path($ref_dir);
my $ref_fa = file($ref_dir, 'pombe_ref.fa')->stringify;
copy($ref_fa_source, $ref_fa);

VRPipe::PipelineSetup->create(
    name       => 'fermikit unitig assembly test',
    datasource => VRPipe::DataSource->create(
        type    => 'fofn_with_metadata',
        method  => 'all',
        source  => file(qw(t data pombe_bam.fofnwm))->absolute->stringify,
        options => {},
    ),
    output_root => $output_dir,
    pipeline    => $assembly_pipeline,
    options     => {
        samtools_exe                            => 'samtools',
        samtools_bam2fq_options                 => '',
        bfc_exe                                 => 'bfc',
        bfc_error_correction_options            => '-s 12m -t 4',
        bfc_filter_options                      => '-1s 12m -k 61 -t 4',
        ropebwt2_exe                            => 'ropebwt2',
        ropebwt2_options                        => '-dNCr',
        fermi2_exe                              => 'fermi2',
        fermi2_assemble_options                 => '-l 71 -m 113 -t 4',
        fermi2_simplify_options                 => '',
        cleanup                                 => 0,
        bfc_error_correct_memory                => 2000,
        bfc_filter_error_corrected_reads_memory => 2000,
        fmd_index_memory                        => 2000,
        fermi2_assemble_memory                  => 2000,
        fermi2_simplify_memory                  => 2000,
    }
);

VRPipe::PipelineSetup->create(
    name       => 'fermikit calling test',
    datasource => VRPipe::DataSource->create(
        type    => 'vrpipe',
        method  => 'all',
        source  => '1[5]',
        options => {}
    ),
    output_root => $output_dir,
    pipeline    => $calling_pipeline,
    options     => {
        reference_fasta                   => $ref_fa,
        reference_assembly_name           => 'pombe1',
        reference_public_url              => 'ftp://s.pombe.com/ref.fa',
        reference_species                 => 'S.pombe',
        bwa_exe                           => 'bwa',
        bwa_index_options                 => '-a is',
        bwa_mem_options                   => '-x intractg',
        htsbox_exe                        => 'htsbox',
        htsbox_pileup_options             => '-Ccu',
        htsbox_abreak_options             => '-bcu',
        fasta_index_memory                => 500,
        sequence_dictionary_memory        => 500,
        bwa_index_memory                  => 500,
        bwa_mem_align_unitigs_memory      => 500,
        bam_sort_memory                   => 500,
        htsbox_pileup_memory              => 500,
        htsbox_abreak_memory              => 500,
        hapdip_filter_pileup_calls_memory => 500,
    }
);

ok handle_pipeline(), 'fermikit pipelines ran ok';

my @assembly_files;
my @output_subdirs = output_subdirs(1, 1);
push(@assembly_files, file(@output_subdirs, '1_bfc_error_correct',                qq[pombe.fq.gz]));
push(@assembly_files, file(@output_subdirs, '2_bfc_filter_error_corrected_reads', qq[pombe.fq.gz]));
push(@assembly_files, file(@output_subdirs, '3_fmd_index',                        qq[pombe.fmd]));
push(@assembly_files, file(@output_subdirs, '4_fermi2_assemble',                  qq[pombe.fq.gz]));
push(@assembly_files, file(@output_subdirs, '5_fermi2_simplify',                  qq[pombe.fq.gz]));

my @ref_files;
foreach my $suffix (qw(fai dict bwt sa amb ann pac)) {
    push @ref_files, file($output_dir, 'ref', 'pombe_ref.fa.' . $suffix);
}

my @calling_files;
@output_subdirs = output_subdirs(2, 2);
push(@calling_files, file(@output_subdirs, '4_bwa_mem_align_unitigs', qq[0.unitig.bam]));
push(@calling_files, file(@output_subdirs, '5_bam_sort',              qq[0.unitig.sorted.bam]));
push(@calling_files, file(@output_subdirs, '6_unitig_pileup_calling', qq[pileup.vcf.gz]));
push(@calling_files, file(@output_subdirs, '6_unitig_pileup_calling', qq[pileup.vcf.gz.tbi]));
push(@calling_files, file(@output_subdirs, '7_unitig_abreak_calling', qq[0.unitig.vcf.gz]));
push(@calling_files, file(@output_subdirs, '7_unitig_abreak_calling', qq[0.unitig.vcf.gz.tbi]));
# push(@calling_files, file(@output_subdirs, '8_filter_unitig_pileup_calls', qq[pileup.vcf.gz]));
# push(@calling_files, file(@output_subdirs, '8_filter_unitig_pileup_calls', qq[pileup.vcf.gz.tbi]));

ok handle_pipeline(@assembly_files, @ref_files, @calling_files), 'fermikit pipelines created expected output files';

finish;
