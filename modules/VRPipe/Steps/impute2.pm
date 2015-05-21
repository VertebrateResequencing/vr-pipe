
=head1 NAME

VRPipe::Steps::impute2 - a step

=head1 DESCRIPTION

This step runs IMPUTE2 on chunks of variants extracted from VCF files given a
set of bed files defining the coordinates of chunks. It converts and splits 
the input VCF on the fly into formats required by impute2. It then runs 
imputation on each chunk and converts the output chunks back to VCF.

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

class VRPipe::Steps::impute2 with VRPipe::StepRole {
    method options_definition {
        return {
            impute2_exe => VRPipe::StepOption->create(
                description   => 'path to your impute2 exe',
                optional      => 1,
                default_value => 'impute2'
            ),
            impute2_options => VRPipe::StepOption->create(
                description   => 'options to impute2 excluding --buffer, -m, -h, -l, -g, -int, -o and -sample_g',
                optional      => 1,
                default_value => '-Ne 20000 -k 80 -allow_large_regions'
            ),
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools exe',
                optional      => 1,
                default_value => 'bcftools'
            ),
            vcf_GL_tag => VRPipe::StepOption->create(
                description   => 'genotype likelihood tag in the VCF (e.g. PL, GL or GT) to be used for imputation',
                optional      => 1,
                default_value => 'GT'
            ),
            tabix_exe => VRPipe::StepOption->create(
                description   => 'path to your tabix exe',
                optional      => 1,
                default_value => 'tabix'
            ),
            ref_vcf => VRPipe::StepOption->create(
                description => 'path to a reference VCF that will be passed to impute2 via -h and -l options (may contain the string "{CHROM}"). When multiple reference panels are given, the pipeline will run impute2 with the -merge_ref_panels option.',
                optional    => 1,
            ),
            genetic_map => VRPipe::StepOption->create(
                description   => 'path to your genetic map files (may contain the string "{CHROM}" which will be expanded)',
                optional      => 0,
                default_value => '/nfs/users/nfs_p/pd3/sandbox/svn/impute2/ALL_1000G_phase1interim_jun2011_impute/genetic_map_chr{CHROM}_combined_b37.txt'
            ),
            bgzip_exe => VRPipe::StepOption->create(
                description   => 'path to the bgzip executable',
                optional      => 1,
                default_value => 'bgzip'
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
            
            my $impute2_exe  = $options->{impute2_exe};
            my $impute2_opts = $options->{impute2_options};
            my $tabix_exe    = $options->{tabix_exe};
            my $bgzip_exe    = $options->{bgzip_exe};
            my $bcftools_exe = $options->{bcftools_exe};
            my $ref_vcf      = $options->{ref_vcf};
            my $vcf_GL_tag   = $options->{vcf_GL_tag};
            
            if ($impute2_opts =~ /\s-(buffer|m|int|h|l|g|sample_g|o)\s/) {
                $self->throw("impute2_options should not include --buffer, -m, -h, -l, -g, -int, -o or -sample_g");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'impute2',
                    version => VRPipe::StepCmdSummary->determine_version($impute2_exe, '^.*version (.+)$'),
                    summary => "impute2 $impute2_opts -buffer 0 -include_buffer_in_output -o_gz -m \$genetic_map_file -int \$start \$end -h \$chunk_haps.gz -l \$chunk_legend.gz -g \$chunk_gts.gz -o \$chunk_output"
                )
            );
            
            my @ref_vcfs = split(/\s+/, $ref_vcf);
            
            my $in_vcf  = @{ $self->inputs->{vcf_files} }[0];
            my $meta    = $in_vcf->metadata;
            my $in_path = $in_vcf->path;
            
            my $header_file = $self->output_file(output_key => 'info_files', basename => "scores.hdr", type => 'txt');
            my $hdr = $header_file->openw;
            print $hdr "##INFO=<ID=IMP2,Number=3,Type=Float,Description=\"IMPUTE2 scores: exp_freq_a1, info, certainty\">";
            $hdr->close;
            my $header_path = $header_file->path;
            
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
                    
                    my @outfiles = ();
                    
                    $self->output_file(sub_dir => $sub_dir, basename => "01.vcf_to_impute2.$from-$to.gen.gz",     type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "01.vcf_to_impute2.$from-$to.samples",    type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-$to.gz",                type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-${to}_samples",         type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-$to.bcf",               type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-$to.bcf.csi",           type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-$to.concat.bcf",        type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-$to.concat.bcf.csi",    type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-${to}_info.bed",        type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-${to}_info.bed.gz",     type => 'bin', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, basename => "02.impute2.$from-${to}_info.bed.gz.tbi", type => 'bin', temporary => 1);
                    
                    my $log_file         = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.impute2.$from-$to.log",              type => 'txt');
                    my $warnings_file    = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.impute2.$from-${to}_warnings",       type => 'txt');
                    my $summary_file     = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.impute2.$from-${to}_summary",        type => 'txt');
                    my $info_file        = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.impute2.$from-${to}_info",           type => 'txt');
                    my $sample_info_file = $self->output_file(sub_dir => $sub_dir, output_key => 'info_files', basename => "02.impute2.$from-${to}_info_by_sample", type => 'txt');
                    
                    my $out_vcf = $self->output_file(sub_dir => $sub_dir, output_key => 'vcf_files', basename => "03.impute2_to_vcf.$from-$to.vcf.gz",     type => 'vcf', metadata => { %$meta, chrom => $chr, from => $from, to => $to });
                    my $out_tbi = $self->output_file(sub_dir => $sub_dir, output_key => 'tbi_files', basename => "03.impute2_to_vcf.$from-$to.vcf.gz.tbi", type => 'bin', metadata => { %$meta, chrom => $chr, from => $from, to => $to });
                    
                    my $out_path = $out_vcf->path;
                    
                    my $gen_map = $self->expand_chrom($options->{genetic_map}, $chr);
                    my $ref_path = $self->expand_chrom($ref_vcf, $chr);
                    
                    push(@outfiles, ($log_file, $warnings_file, $summary_file, $info_file, $sample_info_file, $out_vcf, $out_tbi));
                    
                    ## Time was observed to scale linearly with number of markers (M*2 => T*2) and number of template haplotypes (K*2 => T*2), even though relationship with K is known to be superlinear or quadratic for large K. We estimate the time given it takes 500 sd to run with M=5000 and K=80
                    my $Markers = $chunk[3];
                    my ($N_haps) = ($impute2_opts =~ m/-k\s+/) ? $impute2_opts =~ m/-k\s+(\d+)/ : 80;
                    my $time = int(500 * (1 + ($N_haps - 80) / 80) * (1 + ($Markers - 5000) / 5000)); #seconds
                    
                    my $req = $self->new_requirements(memory => 1000, time => $time);
                    my $this_cmd = "use VRPipe::Steps::impute2; VRPipe::Steps::impute2->run_impute2(chunk => [qw(@chunk)], in_vcf => q[$in_path], out_vcf => q[$out_path], ref_vcf => q[$ref_path], impute2 => q[$impute2_exe], impute2_opts => q[$impute2_opts], bcftools => q[$bcftools_exe], tabix => q[$tabix_exe], bgzip => q[$bgzip_exe], gen_map => q[$gen_map], vcf_GL_tag => q[$vcf_GL_tag], header_file => q[$header_path]);";
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
        return "Run IMPUTE2 to impute or phase chunks of input VCF files";
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
    
    method run_impute2 (ClassName|Object $self: Str|File :$in_vcf!, ArrayRef[Str] :$chunk!, Str|File :$out_vcf!, Str|File :$ref_vcf!, Str :$impute2!, Str :$impute2_opts!, Str :$bcftools!, Str :$tabix!, Str :$bgzip!, Str :$gen_map!, Str :$vcf_GL_tag!, Str :$header_file!) {
        my @region = @{$chunk};
        my $chr    = $region[0];
        my $from   = $region[1];
        my $to     = $region[2];
        
        my $out_file = VRPipe::File->get(path => $out_vcf);
        my $outdir = $out_file->dir;
        
        my $ref   = '';
        my $known = '';
        if ($ref_vcf) {
            my @ref_vcfs = split(/\s+/, $ref_vcf);
            foreach my $ref_path (@ref_vcfs) {
                my $ref_file = VRPipe::File->get(path => $ref_path);
                my $refdir = $ref_file->dir;
                $ref .= " -h $refdir/$chr/$from-${to}.hap.gz -l $refdir/$chr/$from-${to}.legend.gz";
                $self->throw("files $refdir/$chr/$from-${to}* missing") unless -e "$refdir/$chr/$from-$to.ref.sites";
            }
            $ref .= " -merge_ref_panels" unless @ref_vcfs == 1;
        }
        
        my $cmd_line = "$bcftools convert --tag $vcf_GL_tag -r $chr:$from-$to --gensample $outdir/01.vcf_to_impute2.$from-$to $in_vcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $out         = "$outdir/02.impute2.$from-$to";
        my $impute2_cmd = "$impute2 $impute2_opts -buffer 0 -include_buffer_in_output -o_gz -m $gen_map -int $from $to $ref -g $outdir/01.vcf_to_impute2.$from-$to.gen.gz -sample_g $outdir/01.vcf_to_impute2.$from-$to.samples -o $out > $out.log";
        $out_file->disconnect;
        system($impute2_cmd) && $self->throw("failed to run [$impute2_cmd]");
        
        # Sanity check outputs
        my $warnings_file = VRPipe::File->get(path => "${out}_warnings");
        my $info_file     = VRPipe::File->get(path => "${out}_info");
        my $fh            = $warnings_file->openr;
        my @warns         = <$fh>;
        $fh->close;
        $fh = $info_file->openr;
        my @info = <$fh>;
        $fh->close;
        if (@warns > @info * 0.1) { $self->throw("Too many warnings, better to check this: ${out}_warnings"); }
        
        # convert impute2 genotypes to bcf format
        $cmd_line = "$bcftools convert --gensample2vcf $out.gz,${out}_samples -Ob -o $out.bcf && $bcftools index $out.bcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $cmd_line = "$bcftools concat -r $chr:$from-$to -aD $out.bcf $in_vcf -Ob -o $out.concat.bcf && $bcftools index $out.concat.bcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        # convert impute2 info scores to bed format
        $fh = $info_file->openr;
        my $info_bed_file = VRPipe::File->get(path => "${out}_info.bed");
        my $bfh = $info_bed_file->openw;
        while (<$fh>) {
            next if ($_ =~ /certainty/);
            my @scores = split(/\s+/, $_);
            my $str    = $scores[2] - 1;
            my $end    = $scores[2];
            print $bfh "$chr\t$str\t$end\t$scores[3],$scores[4],$scores[5]\n";
        }
        $bfh->close;
        $fh->close;
        $cmd_line = "cat ${out}_info.bed | $bgzip -c > ${out}_info.bed.gz && $tabix -p bed ${out}_info.bed.gz";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $cmd_line = "$bcftools annotate -r $chr:$from-$to -c ID,QUAL,FILTER,INFO,^FMT/GT,^FMT/GP -a $in_vcf -Ou $out.concat.bcf | $bcftools annotate -c CHROM,FROM,TO,IMP2 -h $header_file -a ${out}_info.bed.gz -Oz -o $out_vcf && $tabix -p vcf $out_vcf";
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        return 1;
    }

}

1;
