
=head1 NAME

VRPipe::Steps::sequenom_csv_to_vcf - a step

=head1 DESCRIPTION

This step is Sanger-specific since it uses sequenom plate plex metadata to find
a file in Sanger's filesystem that contains coordinate and strand information
on the SNPs.

It converts the data in a sequenom CSV file into an accurate VCF with correct
SNP positions and strand orientation.

It also considers the gender markers and assigns sequenom_gender metadata.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::sequenom_csv_to_vcf extends VRPipe::Steps::irods {
    use DateTime;
    
    method _build_irods_exes {
        return { iget => 'iget', iquest => 'iquest', ichksum => 'ichksum' };
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            irods_get_zone => VRPipe::StepOption->create(
                description   => 'the zone (top level directory) where your data is stored in iRODs',
                optional      => 1,
                default_value => 'seq'
            ),
            snp_manifest => VRPipe::StepOption->create(
                description => 'file with position and strand of each SNP; if not supplied, an appropriate one will be found and used',
                optional    => 1
            ),
            vcf_sort_exe => VRPipe::StepOption->create(
                description   => 'path to the vcf-sort executable',
                optional      => 1,
                default_value => 'vcf-sort'
            ),
            vcf_sort_options => VRPipe::StepOption->create(
                description   => 'options to vcf-sort',
                optional      => 1,
                default_value => '--chromosomal-order'
            ),
            bgzip_exe => VRPipe::StepOption->create(
                description   => 'path to the bgzip executable',
                optional      => 1,
                default_value => 'bgzip'
            )
        };
    }
    
    method inputs_definition {
        return {
            csv_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'csv file with sequenom_plex metadata',
                max_files   => -1,
                # metadata    => { sequenom_plex => 'sequenom plate plex identifier' }
            ),
        };
    }
    
    method outputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'VCF file',
                max_files   => -1,
                metadata    => { sequenom_gender => 'gender of sample derived from sequenom result' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $vcf_sort      = $options->{vcf_sort_exe};
            my $vcf_sort_opts = $options->{vcf_sort_options};
            if ($vcf_sort_opts) {
                $vcf_sort .= ' ' . $vcf_sort_opts;
            }
            my $bgzip = $options->{bgzip_exe};
            my $snp_manifest_file = $options->{snp_manifest} || '';
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $csv_file (@{ $self->inputs->{csv_files} }) {
                my $basename = $csv_file->basename;
                $basename =~ s/\.csv$//;
                my $unsorted = $basename . '.unsorted.vcf';
                $basename .= '.vcf.gz';
                $self->output_file(basename => $unsorted, type => 'vcf', temporary => 1);
                my $vcf_file_path = $self->output_file(output_key => 'vcf_files', basename => $basename, type => 'vcf', metadata => $csv_file->metadata)->path;
                my $csv_file_path = $csv_file->path;
                
                my $cmd = "use VRPipe::Steps::sequenom_csv_to_vcf; VRPipe::Steps::sequenom_csv_to_vcf->csv_to_vcf(csv => q[$csv_file_path], vcf => q[$vcf_file_path], snp_manifest => q[$snp_manifest_file], vcf_sort => q[$vcf_sort], bgzip => q[$bgzip]);";
                $self->dispatch_vrpipecode($cmd, $req);
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
        return "Convert a CSV file containing sequenom results into a sorted compressed VCF suitable for calling with.";
    }
    
    method csv_to_vcf (ClassName|Object $self: Str|File :$csv, Str|File :$vcf, Str :$vcf_sort, Str :$bgzip, Maybe[Str|File] :$snp_manifest?) {
        my $csv_file = VRPipe::File->get(path => $csv);
        my $vcf_file = VRPipe::File->get(path => $vcf);
        my $vcf_file_unsorted = $vcf;
        $vcf_file_unsorted =~ s/\.vcf.gz$/.unsorted.vcf/;
        $vcf_file_unsorted = VRPipe::File->get(path => $vcf_file_unsorted);
        
        unless ($snp_manifest) {
            # search for it in irods based on the $csv sequenom_plex
            my $sequenom_plex = $csv_file->meta_value('sequenom_plex');
            #***...
            die "determining snp manifest from sequenom_plex not yet implemented\n";
        }
        
        open(my $mfh, '<', $snp_manifest) || die "Could not read from $snp_manifest\n";
        my $ofh = $vcf_file_unsorted->openw;
        
        # get date and sample name for VCF header
        my $dt     = DateTime->now;
        my $date   = $dt->ymd('');
        my $sample = $csv_file->meta_value('sample');
        $vcf_file->disconnect;
        
        # print VCF header
        print $ofh "##fileformat=VCFv4.0\n";
        print $ofh "##fileDate=$date\n";
        print $ofh "##source=$csv sequenom results\n";
        print $ofh "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
        print $ofh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        print $ofh "##FORMAT=<ID=IA,Number=1,Type=Float,Description=\"Intensity of the A Allele\">\n";
        print $ofh "##FORMAT=<ID=IB,Number=1,Type=Float,Description=\"Intensity of the B Allele\">\n";
        print $ofh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";
        
        # read the snp manifest data into a hash
        my %manifest;
        <$mfh>; # header
        while (<$mfh>) {
            chomp;
            # SNP_NAME REF_ALLELE ALT_ALLELE CHR POS SEQUENOM_ASSAY_STRAND
            my ($snp_name, @vals) = split(/\t/);
            
            if (defined $manifest{$snp_name}) {
                # SNP_NAME can be repeated twice in the case that the SNP is on
                # both X and Y because we're dealing with a gender marker
                $manifest{$snp_name}->[1] = $vals[0];
                push(@{ $manifest{$snp_name} }, 1);
            }
            else {
                $manifest{$snp_name} = \@vals;
            }
        }
        close($mfh);
        
        # go through the csv file and output the VCF records
        my $ifh = $csv_file->openr;
        $vcf_file->disconnect;
        my %gender;
        my $records = 0;
        <$ifh>; # header
        while (<$ifh>) {
            my ($allele, $assay_id, undef, undef, undef, $genotype, $height) = split(/\t/);
            my $second_line = <$ifh>;
            $second_line || die "Uneven number of lines; no partner for $assay_id found\n";
            my ($allele_two, $assay_id_two, undef, undef, undef, undef, $height_two) = split(/\t/, $second_line);
            die "Mismatching lines: $assay_id_two followed $assay_id\n" unless $assay_id eq $assay_id_two;
            $assay_id =~ s/^W\d+\-//;
            $height =~ s/^\./0./;
            $height_two =~ s/^\./0./;
            
            if (length($genotype) == 1) {
                $allele     = $genotype;
                $allele_two = $genotype;
            }
            
            # get the snp meta info from the manifest hash
            my ($ref_allele, $alt_allele, $chr, $pos, $strand, $gender_marker) = @{ $manifest{$assay_id} || next };
            
            # flip to the SEQUENOM_ASSAY_STRAND strand if necessary
            if ($strand eq '-') {
                $ref_allele =~ tr/ATGCatgc/TACGtacg/;
                $alt_allele =~ tr/ATGCatgc/TACGtacg/;
            }
            
            # match the alleles against the ref alleles
            my $match_one = my $match_two = 0;
            if ($allele ne $ref_allele) {
                $match_one = 1;
            }
            if ($allele_two ne $ref_allele) {
                $match_two = 1;
            }
            
            #*** no idea if further flipping as per genome_studio_fcr_to_vcf
            #    step is required
            
            # if this was for a gender marker, record the result
            if ($gender_marker) {
                if (length($genotype) == 2) {
                    $gender{M}++;
                }
                else {
                    $gender{F}++;
                }
            }
            
            # write conversion
            print $ofh "$chr\t$pos\t$assay_id\t$ref_allele\t$alt_allele\t.\t.\tNS=1\tGT:IA:IB\t$match_one/$match_two:$height:$height_two\n";
            $records++;
        }
        $csv_file->close;
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
        
        # record the gender, defaulting to male if in doubt
        my $gender;
        my @genders = keys %gender;
        if (@genders) {
            if (@genders == 1) {
                $gender = $genders[0];
            }
            else {
                $gender = 'M';
            }
        }
        else {
            $gender = 'M';
        }
        $vcf_file->add_metadata({ sequenom_gender => $gender });
    }
}

1;
