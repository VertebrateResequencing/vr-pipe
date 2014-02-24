
=head1 NAME

VRPipe::Steps::genome_studio_fcr_to_vcf - a step

=head1 DESCRIPTION

This step converts the data in a single-sample fcr file (eg. produced by the
split_genome_studio_genotype_files step) into an accurate VCF with correct SNP
positions and strand orientation.

It uses the fcr-to-vcf exe from the vr-codebase git repository.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2014 Genome Research Limited.

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

class VRPipe::Steps::genome_studio_fcr_to_vcf with VRPipe::StepRole {
    use DateTime;
    
    method options_definition {
        return {
            vcf_sample_from_metadata => VRPipe::StepOption->create(
                description   => 'if the id in the second column of the fcr matches metadata with key x, but you want the sample id in the VCF to have metadata from key y, provide x:y',
                optional      => 1,
                default_value => 'library:sample'
            ),
            fcr_to_vcf_exe => VRPipe::StepOption->create(
                description   => 'path to the fcr-to-vcf executable',
                optional      => 1,
                default_value => 'fcr-to-vcf'
            )
        };
    }
    
    method inputs_definition {
        return {
            map_file => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'mapping file produced by the illumina_coreexome_manifest_to_map step'
            ),
            fcr_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'single-sample fcr file with sample metadata',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line' }
            ),
        };
    }
    
    method outputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'VCF file',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self       = shift;
            my $options    = $self->options;
            my $fcr_to_vcf = $options->{fcr_to_vcf_exe};
            my $sfm        = $options->{vcf_sample_from_metadata};
            my ($src_key, $dst_key);
            my $summary_s = '';
            if ($sfm) {
                ($src_key, $dst_key) = split(':', $sfm);
                $summary_s = " -s $sfm";
            }
            
            my ($map_file) = @{ $self->inputs->{map_file} };
            $map_file = $map_file->path;
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'fcr-to-vcf',
                    version => 0,
                    summary => "cat \$fcr_file | $fcr_to_vcf -a \$map_file$summary_s -o \$outdir"
                )
            );
            
            foreach my $fcr_file (@{ $self->inputs->{fcr_files} }) {
                my $meta = $fcr_file->metadata;
                
                my $outdir = 'output';
                my $s      = '';
                if ($src_key && $dst_key && defined $meta->{$src_key} && defined $meta->{$dst_key}) {
                    $s      = " -s $meta->{$src_key}:$meta->{$dst_key}";
                    $outdir = $meta->{$dst_key};
                }
                
                #$self->output_file(basename => $basename, type => 'vcf', temporary => 1);
                
                my $vcf_file_path = $self->output_file(output_key => 'vcf_files', basename => file($outdir, "$outdir.vcf.gz"), type => 'vcf', metadata => $meta)->path;
                my $fcr_file_path = $fcr_file->path;
                
                $self->dispatch(["cat $fcr_file_path | $fcr_to_vcf -a $map_file$s -o $outdir", $req]);
            }
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method description {
        return "Convert single-sample genotype data text files (FCR files) to sorted compressed VCF suitable for calling with.";
    }
    
    method fcr_to_vcf (ClassName|Object $self: Str|File :$fcr, Str|File :$vcf, Str|File :$snp_manifest, Str :$vcf_sort, Str :$bgzip) {
        open(my $mfh, '<', $snp_manifest) || die "Could not read from $snp_manifest\n";
        my $fcr_file = VRPipe::File->get(path => $fcr);
        my $vcf_file = VRPipe::File->get(path => $vcf);
        my $vcf_file_unsorted = $vcf;
        $vcf_file_unsorted =~ s/\.vcf.gz$/.unsorted.vcf/;
        $vcf_file_unsorted = VRPipe::File->get(path => $vcf_file_unsorted);
        my $ofh = $vcf_file_unsorted->openw;
        
        # get date and sample name for VCF header; the fcr file may have the
        # library name instead of sample name, so get that as well
        my $dt      = DateTime->now;
        my $date    = $dt->ymd('');
        my $sample  = $fcr_file->meta_value('sample');
        my $library = $fcr_file->meta_value('library');
        $vcf_file->disconnect;
        
        # print VCF header
        print $ofh "##fileformat=VCFv4.0\n";
        print $ofh "##fileDate=$date\n";
        print $ofh "##source=$fcr HipSci genotyping file\n";
        print $ofh "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
        print $ofh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        print $ofh "##FORMAT=<ID=GC,Number=1,Type=Float,Description=\"GenCall score\">\n";
        print $ofh "##FORMAT=<ID=IA,Number=1,Type=Float,Description=\"Intensity of the A Allele\">\n";
        print $ofh "##FORMAT=<ID=IB,Number=1,Type=Float,Description=\"Intensity of the B Allele\">\n";
        print $ofh "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">\n";
        print $ofh "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">\n";
        print $ofh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";
        
        # read the snp manifest data into a hash
        my %manifest;
        while (<$mfh>) {
            chomp;
            my @vals = split(/\t/);
            next if $vals[5] eq 'UNKNOWN'; # (this happens if the strand file did not have entries for everything in raw manifest)
            
            my $manifestSummarySNP = $vals[2];
            $manifestSummarySNP =~ /\[([A-Za-z]+)\/([A-Za-z]+)/;
            my $firstRefAllele  = $1;
            my $secondRefAllele = $2;
            
            $manifest{ $vals[0] } = [$vals[1], $firstRefAllele, $secondRefAllele, $vals[3], $vals[4], $vals[5]];
        }
        close($mfh);
        
        # go through the fcr file and output the VCF records
        my $ifh = $fcr_file->openr;
        $vcf_file->disconnect;
        my $records = 0;
        while (<$ifh>) {
            chomp;
            my @vals = split(/\t/);
            next unless $vals[1] eq $library || $vals[1] eq $sample;
            
            my $genotypeSNPName   = $vals[0];
            my $genotypeAlleleOne = $vals[2];
            my $genotypeAlleleTwo = $vals[3];
            my $genotypeGCScore   = $vals[4];
            
            # skip non-SNPs
            if (($genotypeAlleleOne eq 'D') || ($genotypeAlleleOne eq 'I')) {
                next;
            }
            if (($genotypeAlleleOne eq '-') || ($genotypeAlleleTwo eq '-')) {
                next;
            }
            
            my $intensityA = $vals[7];
            my $intensityB = $vals[8];
            my $BAF        = $vals[11];
            my $LRR        = $vals[12];
            
            # get the snp meta info from the manifest hash
            my ($manifestSummaryIlmnStrand, $firstRefAllele, $secondRefAllele, $manifestSummaryChr, $manifestSummaryMapInfo, $manifestSummarySourceStrand) = @{ $manifest{$genotypeSNPName} || next };
            
            # flip to the ILMN strand if necessary (ie. 'BOT' in manifest)
            if ($manifestSummaryIlmnStrand eq 'BOT') {
                # flip to complements
                $firstRefAllele =~ tr/ATGCatgc/TACGtacg/;
                $secondRefAllele =~ tr/ATGCatgc/TACGtacg/;
            }
            
            # match the alleles against the ref alleles
            my $matchOne = my $matchTwo = 0;
            if ($genotypeAlleleOne ne $firstRefAllele) {
                $matchOne = 1;
            }
            if ($genotypeAlleleTwo ne $firstRefAllele) {
                $matchTwo = 1;
            }
            
            # flip to the forward strand if necessary (ie. '-' in manifest)
            if ($manifestSummarySourceStrand eq '-') {
                $firstRefAllele =~ tr/ATGCatgc/TACGtacg/;
                $secondRefAllele =~ tr/ATGCatgc/TACGtacg/;
            }
            
            # write conversion
            print $ofh "$manifestSummaryChr\t$manifestSummaryMapInfo\t$genotypeSNPName\t$firstRefAllele\t$secondRefAllele\t.\t.\tNS=1\tGT:GC:IA:IB:BAF:LRR\t$matchOne/$matchTwo:$genotypeGCScore:$intensityA:$intensityB:$BAF:$LRR\n";
            $records++;
        }
        $fcr_file->close;
        $vcf_file_unsorted->close;
        
        # sort the vcf file
        my $vcf_sort_cmd = "$vcf_sort " . $vcf_file_unsorted->path . " | $bgzip -c > " . $vcf_file->path;
        system($vcf_sort_cmd) && die "VCF sort [$vcf_sort_cmd] failed\n";
        
        # check it has the correct number of lines
        $vcf_file->update_stats_from_disc;
        my $actual_records = $vcf_file->num_records;
        unless ($actual_records == $records) {
            die "Expected $records records in the sorted VCF file, but only got $actual_records\n";
        }
    }
}

1;
