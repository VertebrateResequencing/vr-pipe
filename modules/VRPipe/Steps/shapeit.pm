
=head1 NAME

VRPipe::Steps::shapeit - a step

=head1 DESCRIPTION

This step runs SHAPEIT on chunks of variants extracted from VCF files given a
set of bed files defining the coordinates of chunks. It converts and splits 
the input VCF on the fly into formats required by shapeit. It then runs the
phasing algorithm on each chunk and converts the shapeit output back to VCF. 
The step's output will be phased chunks in VCF format.

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk> and Petr Danecek <pd3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Steps::shapeit with VRPipe::StepRole {
    method options_definition {
        return {
            shapeit_exe => VRPipe::StepOption->create(
                description   => 'path to your shapeit exe',
                optional      => 1,
                default_value => 'shapeit'
            ),
            shapeit_options => VRPipe::StepOption->create(
                description => 'options to shapeit, exclude options that concern input/output file definition',
                optional    => 1,
            ),
            tabix_exe => VRPipe::StepOption->create(
                description   => 'path to your tabix exe',
                optional      => 1,
                default_value => 'tabix'
            ),
            ref_vcf => VRPipe::StepOption->create(
                description => 'path to a reference VCF which may contain the string "{CHROM}" that will be expanded to chromosome names (if multiple VCFs are specified, only the first one is used with shapeit)',
                optional    => 1,
            ),
            genetic_map => VRPipe::StepOption->create(
                description   => 'path to your genetic map files (may contain the string "{CHROM}" which will be expanded)',
                optional      => 0,
                default_value => '/nfs/users/nfs_p/pd3/sandbox/svn/impute2/ALL_1000G_phase1interim_jun2011_impute/genetic_map_chr{CHROM}_combined_b37.txt'
            ),
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools exe',
                optional      => 1,
                default_value => 'bcftools'
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                min_files   => 1,
                max_files   => 1,
                description => '1 VCF file'
            ),
            bed_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => -1,
                description => '1 or more bed files'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $shapeit_exe  = $options->{shapeit_exe};
            my $shapeit_opts = $options->{shapeit_options};
            my $tabix_exe    = $options->{tabix_exe};
            my $bcftools_exe = $options->{bcftools_exe};
            
            if ($shapeit_opts =~ /\s--(input-|output-|include-|thread)\s/) {
                $self->throw("shapeit_options should not include --input-* or --output-* or --include-* or --thread");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'shapeit',
                    version => VRPipe::StepCmdSummary->determine_version($shapeit_exe, 'Version : (.+)$'),
                    summary => "shapeit $shapeit_opts --input-ref \$chunk_haps.gz \$chunk_legend.gz \$chunk_samples --include-snp \$chunk.ref.sites --input-map \$genetic_map_file --input-gen \$chunk_haps.gz \$chunk_samples --output-max \$chunk.haps.gz \$chunk.samples --output-log \$chunk.log"
                )
            );
            
            my $in_vcf  = @{ $self->inputs->{vcf_files} }[0];
            my $in_path = $in_vcf->path;
            my @samples = grep { chomp } `$bcftools_exe query -l $in_path`;
            
            my $meta = $in_vcf->metadata;
            
            foreach my $bed_file (@{ $self->inputs->{bed_files} }) {
                my $fh = $bed_file->openr;
                while (<$fh>) {
                    chomp;
                    next if ($_ =~ /^#/);
                    my @chunk   = split(/\t/, $_);
                    my $chr     = $chunk[0];
                    my $from    = $chunk[1];
                    my $to      = $chunk[2];
                    my $sub_dir = $chr;
                    
                    $self->output_file(sub_dir => $sub_dir, basename => "01.vcf_to_shapeit.$from-$to.gen.gz",  type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "01.vcf_to_shapeit.$from-$to.samples", type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.shapeit.$from-$to.haps.gz",        type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.shapeit.$from-$to.samples",        type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.shapeit.$from-$to.bcf",            type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.shapeit.$from-$to.bcf.csi",        type => 'bin', temporary => 1);
                    
                    my $log_file       = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.shapeit.$from-$to.log",                type => 'txt');
                    my $snp_mm_file    = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.shapeit.$from-$to.snp.mm",             type => 'txt');
                    my $ind_mm_file    = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.shapeit.$from-$to.ind.mm",             type => 'txt');
                    my $strand_file    = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.shapeit.$from-$to.snp.strand",         type => 'txt');
                    my $exclusion_file = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.shapeit.$from-$to.snp.strand.exclude", type => 'txt');
                    
                    my $out_vcf = $self->output_file(sub_dir => $sub_dir, output_key => 'vcf_files', basename => "03.shapeit_to_vcf.$from-$to.vcf.gz",     type => 'vcf', metadata => { %$meta, chrom => $chr, from => $from, to => $to });
                    my $out_tbi = $self->output_file(sub_dir => $sub_dir, output_key => 'tbi_files', basename => "03.shapeit_to_vcf.$from-$to.vcf.gz.tbi", type => 'bin', metadata => { %$meta, chrom => $chr, from => $from, to => $to });
                    
                    my $out_path = $out_vcf->path;
                    my $ref_path = $self->expand_chrom($options->{ref_vcf}, $chr);
                    my $gen_map  = $self->expand_chrom($options->{genetic_map}, $chr);
                    
                    ## Shapeit, if run in chunks, is more cpu bound than memory bound! Time is expected to scale linearly, and was tested to scale with number of markers (M*2 => T*1.5), samples size (S*40 => T*160) and number of conditioning haplotypes (N*2 => T*1.5). We estimate the time given it takes 200 sd to run with M=5000 and S=100 and N=100
                    my $Markers = $chunk[3];
                    my ($N_haps) = ($shapeit_opts =~ m/--states/) ? $shapeit_opts =~ m/--states\s*(\d+)/ : 100;
                    my $samples = scalar @samples;
                    my $time = int(200 * (1 + 0.5 * ($N_haps - 100) / 100) * (1 + 4 * ($samples - 100) / 100) * (1 + 0.5 * ($Markers - 5000) / 5000)); #seconds
                    
                    my $req = $self->new_requirements(memory => 1000, time => $time);
                    my @outfiles = ($log_file, $snp_mm_file, $ind_mm_file, $strand_file, $exclusion_file, $out_vcf, $out_tbi);
                    my $this_cmd = "use VRPipe::Steps::shapeit; VRPipe::Steps::shapeit->run_shapeit(chunk => [qw(@chunk)], in_vcf => q[$in_path], out_vcf => q[$out_path], ref_vcf => q[$ref_path], shapeit => q[$shapeit_exe], shapeit_opts => q[$shapeit_opts], bcftools => q[$bcftools_exe], tabix => q[$tabix_exe], gen_map => q[$gen_map]);";
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
                
                }
                close($fh);
            
            }
        };
    }
    
    method outputs_definition {
        return {
            info_files => VRPipe::StepIODefinition->create(
                type            => 'txt',
                min_files       => 0,
                max_files       => -1,
                description     => 'additional info files',
                check_existence => 0
            ),
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                min_files   => 1,
                max_files   => -1,
                description => 'output phased VCF chunks',
            ),
            tbi_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                min_files   => 1,
                max_files   => -1,
                description => 'output VCF index files',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run SHAPEIT to phase or refine genotypes of chunks of input VCF files";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method expand_chrom (Str $path, Str $region) {
        if (!defined $region) { return $path; }
        $region =~ s/:.*$//;
        $path =~ s/{CHROM}/$region/g;
        return $path;
    }
    
    method run_shapeit (ClassName|Object $self: Str|File :$in_vcf!, ArrayRef[Str] :$chunk!, Str|File :$out_vcf!, Str|File :$ref_vcf!, Str :$shapeit!, Str :$shapeit_opts!, Str :$bcftools!, Str :$tabix!, Str :$gen_map!) {
        my @region = @{$chunk};
        my $chr    = $region[0];
        my $from   = $region[1];
        my $to     = $region[2];
        
        my $out_file = VRPipe::File->get(path => $out_vcf);
        my $outdir = $out_file->dir;
        
        my $ref = '';
        if ($ref_vcf) {
            my @ref_vcfs = split(/\s+/, $ref_vcf);
            my $ref_file = VRPipe::File->get(path => $ref_vcfs[0]);
            my $refdir = $ref_file->dir;
            $ref = "--input-ref $refdir/$chr/$from-${to}.hap.gz $refdir/$chr/$from-${to}.legend.gz $refdir/$chr/$from-${to}.samples";
            $ref .= " --include-snp $refdir/$chr/$from-$to.ref.sites";
            $self->throw("files $refdir/$chr/$from-${to}* missing") unless -e "$refdir/$chr/$from-$to.ref.sites";
        }
        
        my $cmd_line = "$bcftools convert -r $chr:$from-$to --gensample $outdir/01.vcf_to_shapeit.$from-$to $in_vcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $out         = "$outdir/02.shapeit.$from-$to";
        my $shapeit_cmd = "$shapeit $shapeit_opts $ref --input-map $gen_map --input-gen $outdir/01.vcf_to_shapeit.$from-$to.gen.gz $outdir/01.vcf_to_shapeit.$from-$to.samples --output-max $out.haps.gz $out.samples --output-log $out.log";
        $out_file->disconnect;
        my $failed = system($shapeit_cmd);
        # if shapeit fails in first attempt due to strand mismatch, then rerun with --exclude-snp option
        if ($failed) {
            if (-e "$out.snp.strand.exclude") {
                $ref .= " --exclude-snp $out.snp.strand.exclude";
                $shapeit_cmd = "$shapeit $shapeit_opts $ref --input-map $gen_map --input-gen $outdir/01.vcf_to_shapeit.$from-$to.gen.gz $outdir/01.vcf_to_shapeit.$from-$to.samples --output-max $out.haps.gz $out.samples --output-log $out.log";
                $out_file->disconnect;
                system($shapeit_cmd) && $self->throw("failed to run [$shapeit_cmd]");
            }
            else {
                $self->throw("failed to run [$shapeit_cmd]");
            }
        }
        
        $cmd_line = "$bcftools convert --hapsample2vcf $out.haps.gz,$out.samples -Ob -o $out.bcf && $bcftools index $out.bcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $cmd_line = "$bcftools annotate -r $chr:$from-$to -c -FMT/GT -a $out.bcf $in_vcf -Oz -o $out_vcf && $tabix -p vcf $out_vcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        return 1;
    }

}

1;
