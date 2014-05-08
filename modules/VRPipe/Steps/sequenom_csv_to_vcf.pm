
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

class VRPipe::Steps::sequenom_csv_to_vcf extends VRPipe::Steps::irods {
    use DateTime;
    use VRPipe::Persistent::InMemory;
    
    method _build_irods_exes {
        return { iget => 'iget', ichksum => 'ichksum', imeta => 'imeta' };
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            irods_get_zone => VRPipe::StepOption->create(
                description   => 'the zone (top level directory) where your data is stored in iRODs',
                optional      => 1,
                default_value => 'seq'
            ),
            vcf_sample_from_metadata => VRPipe::StepOption->create(
                description   => 'in the output VCF use the sample name from a metadata value stored on the input csv file; separate multiple keys with + symbols - values will be joined with underscores',
                optional      => 1,
                default_value => 'public_name+sample'
            ),
            sequenom_plex_storage_dir    => VRPipe::StepOption->create(description => 'absolute path to a directory where plex manifest files can be stored'),
            sequencescape_reference_name => VRPipe::StepOption->create(
                description   => 'the name of the reference found in sequencescape that future sequencing data would be mapped with',
                default_value => 'Homo_sapiens (1000Genomes)'
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
            ),
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to the bcftools executable',
                optional      => 1,
                default_value => 'bcftools'
            )
        };
    }
    
    method inputs_definition {
        return {
            csv_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'csv file with sequenom_plex metadata',
                max_files   => -1,
                metadata    => {
                    sequenom_plex => 'sequenom plate plex identifier',
                    optional      => ['sequenom_plex']
                }
            ),
        };
    }
    
    method outputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'VCF file',
                max_files   => -1,
                metadata    => {
                    sequenom_plex   => 'sequenom plate plex identifier',
                    sequenom_gender => 'gender of sample derived from sequenom result'
                }
            ),
            vcf_index => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'index of the vcf file',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $bcftools_exe  = $options->{bcftools_exe};
            my $vcf_sort      = $options->{vcf_sort_exe};
            my $vcf_sort_opts = $options->{vcf_sort_options};
            if ($vcf_sort_opts) {
                $vcf_sort .= ' ' . $vcf_sort_opts;
            }
            my $bgzip                    = $options->{bgzip_exe};
            my $imeta                    = $options->{imeta_exe};
            my $iget                     = $options->{iget_exe};
            my $ichksum                  = $options->{ichksum_exe};
            my $zone                     = $options->{irods_get_zone};
            my $vcf_sample_from_metadata = $options->{vcf_sample_from_metadata};
            
            my $manifest_dir = $options->{sequenom_plex_storage_dir};
            $self->throw("sequenom_plex_storage_path does not exist") unless -d $manifest_dir;
            my $ref_name = $options->{sequencescape_reference_name};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $csv_file (@{ $self->inputs->{csv_files} }) {
                my $basename = $csv_file->basename;
                $basename =~ s/\.csv$//;
                my $unsorted = $basename . '.unsorted.vcf';
                $basename .= '.vcf.gz';
                $self->output_file(basename => $unsorted, type => 'vcf', temporary => 1);
                my $vcf_file_path = $self->output_file(output_key => 'vcf_files', basename => $basename, type => 'vcf', metadata => $csv_file->metadata)->path;
                my $csv_file_path = $csv_file->path;
                $self->output_file(output_key => 'vcf_index', basename => $basename . '.csi', type => 'bin')->path;
                
                my $cmd = "use VRPipe::Steps::sequenom_csv_to_vcf; VRPipe::Steps::sequenom_csv_to_vcf->csv_to_vcf(csv => q[$csv_file_path], vcf => q[$vcf_file_path], manifest_dir => q[$manifest_dir], reference_name => q[$ref_name], vcf_sample_from_metadata => q[$vcf_sample_from_metadata], vcf_sort => q[$vcf_sort], bgzip => q[$bgzip], bcftools => q[$bcftools_exe], imeta => q[$imeta], iget => q[$iget], ichksum => q[$ichksum], zone => q[$zone]);";
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
    
    method csv_to_vcf (ClassName|Object $self: Str|File :$csv, Str|File :$vcf, Str :$vcf_sort, Str :$bgzip, Str :$bcftools, Str|Dir :$manifest_dir, Str :$reference_name, Str :$imeta, Str :$iget, Str :$ichksum, Str :$zone, Str :$vcf_sample_from_metadata = 'public_name+sample') {
        my $csv_file = VRPipe::File->get(path => $csv);
        my $vcf_file = VRPipe::File->get(path => $vcf);
        my $vcf_file_unsorted = $vcf;
        $vcf_file_unsorted =~ s/\.vcf.gz$/.unsorted.vcf/;
        $vcf_file_unsorted = VRPipe::File->get(path => $vcf_file_unsorted);
        my @vsfm_keys = split(/\+/, $vcf_sample_from_metadata);
        
        # figure out the correct snp manifest file to use based on plex and
        # reference
        my $sequenom_plex = $csv_file->meta_value('sequenom_plex');
        unless ($sequenom_plex) {
            # parse the csv to figure out the plex
            my $fh = $csv_file->openr;
            <$fh>;
            my $line = <$fh>;
            my (undef, $id) = split(/\s+/, $line);
            ($sequenom_plex) = split('-', $id);
            $sequenom_plex || $self->throw("Could not parse out sequenom_plex from $id in " . $csv_file->path);
        }
        
        my $reference_base = $reference_name;
        $reference_base =~ s/[\s\(\)]+//g;
        my $snp_manifest_basename = $sequenom_plex . '.' . $reference_base . '.manifest';
        my $snp_manifest = VRPipe::File->create(path => file($manifest_dir, $snp_manifest_basename));
        unless ($snp_manifest->s) {
            # search for it in irods based on the sequenom_plex and reference
            my ($irods_path) = $self->search_by_metadata(metadata => { reference_name => $reference_name, sequenom_plex => $sequenom_plex }, imeta => $imeta, zone => $zone);
            $irods_path || $self->throw("Could not find a manifest file in zone $zone matching metadata reference_name => $reference_base, sequenom_plex => $sequenom_plex");
            
            # we could have multiple processes in parallel all trying to get
            # this same file at the same time; lock and block
            my $im       = VRPipe::Persistent::InMemory->new;
            my $lock_key = 'sequenom_csv_to_vcf.' . $snp_manifest->path;
            $im->block_until_locked($lock_key);
            $snp_manifest->reselect_values_from_db;
            $snp_manifest->update_stats_from_disc;
            unless ($snp_manifest->s) {
                $im->maintain_lock($lock_key);
                warn "getting ", $snp_manifest->path, " at ", time(), "\n";
                $self->get_file(source => $irods_path, dest => $snp_manifest->path, iget => $iget, ichksum => $ichksum);
            }
            $im->unlock($lock_key);
        }
        
        my $mfh = $snp_manifest->openr;
        my $ofh = $vcf_file_unsorted->openw;
        
        # get date and sample name for VCF header
        my $dt     = DateTime->now;
        my $date   = $dt->ymd('');
        my $sample = join('_', map { $csv_file->meta_value($_) } @vsfm_keys);
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
            if (length($genotype) == 0) {
                $match_one = $match_two = '.';
            }
            else {
                if ($allele ne $ref_allele) {
                    $match_one = 1;
                }
                if ($allele_two ne $ref_allele) {
                    $match_two = 1;
                }
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
        $vcf_file->add_metadata({ sequenom_gender => $gender, sequenom_plex => $sequenom_plex });
        
        # index it
        system("$bcftools index $vcf") && die "Failed to index $vcf\n";
    }
}

1;
