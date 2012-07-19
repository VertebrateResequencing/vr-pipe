
=head1 NAME

VRPipe::Steps::bam_stats - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::bam_stats with VRPipe::StepRole {
    use VRPipe::Parser;
    use VRPipe::Parser::fastq;
    use VRPipe::Utils::Math;
    
    method options_definition {
        return { release_date   => VRPipe::StepOption->create(description => 'for DCC-style filenames, provide the release date',                                                                                            optional => 1),
                 sequence_index => VRPipe::StepOption->create(description => 'for DCC-style filenames and using input bams with poor headers, provide a DCC sequence.index',                                                 optional => 1),
                 rg_from_pu     => VRPipe::StepOption->create(description => 'boolean which if true means that bam headers have their RG identifiers in the PU field instead of ID field of the RG lines in the bam header', optional => 1, default_value => 0) };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to calculate stats for') };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options    = $self->options;
            my $extra_args = '';
            foreach my $extra (qw(release_date sequence_index rg_from_pu)) {
                if ($options->{$extra}) {
                    $extra_args .= ", $extra => q[$options->{$extra}]";
                }
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bas = $self->output_file(
                    output_key => 'bas_files',
                    #output_dir => $bam->dir, # *** important to have them in the same dir as input bam?
                    basename => $bam->basename . '.bas',
                    type     => 'txt',
                    metadata => $bam->metadata);
                my $bam_path = $bam->path;
                my $bas_path = $bas->path;
                my $this_cmd = "use VRPipe::Steps::bam_stats; VRPipe::Steps::bam_stats->bas(q[$bam_path], q[$bas_path]$extra_args);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$bas] });
            }
        };
    }
    
    method outputs_definition {
        return { bas_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'a .bas file for each input bam file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Calculates various statistics about bam files, producing .bas files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    our %tech_to_platform = (SLX       => 'ILLUMINA',
                             '454'     => 'LS454',
                             SOLID     => 'ABI_SOLID',
                             ILLUMINA  => 'ILLUMINA',
                             LS454     => 'LS454',
                             ABI_SOLID => 'ABI_SOLID');

=head2 bam_statistics
 
 Title   : bam_statistics
 Usage   : my %rg_stats = VRPipe::Steps::bam_stats->bam_statistics('in.bam');
 Function: Calculate all the stats per-readgroup for a bam needed for the bas
           format.
 Returns : hash with keys as readgroup ids and values as hash refs. Those refs
           have the keys:
           total_bases, mapped_bases, total_reads, mapped_reads,
           mapped_reads_paired_in_seq, mapped_reads_properly_paired,
           percent_mismatch, avg_qual, avg_isize, sd_isize, median_isize, mad,
           duplicate_reads and duplicate_bases
 Args    : input bam file (must have been processed by samtools fillmd so that
           there are accurate NM tags for each record)

=cut
    
    method bam_statistics (ClassName|Object $self: Str|File $bam_file!) {
        # go through the bam and accumulate all the raw stats in little memory
        my $pb = VRPipe::Parser->create('bam', { file => file($bam_file) });
        $pb->get_fields('SEQ_LENGTH', 'MAPPED_SEQ_LENGTH', 'FLAG', 'QUAL', 'MAPQ', 'ISIZE', 'RG', 'NM');
        my $pr = $pb->parsed_record();
        
        my %readgroup_data;
        my $previous_rg = 'unknown_readgroup';
        while ($pb->next_record) {
            my $rg   = $pr->{RG};
            my $flag = $pr->{FLAG};
            
            unless (defined $rg) {
                $self->warn("a read had no RG tag, using previous RG tag '$previous_rg'");
                $rg = $previous_rg;
            }
            $previous_rg = $rg;
            
            my @this_rg_data = @{ $readgroup_data{$rg} || [] };
            $this_rg_data[0]++;
            $this_rg_data[1] += $pr->{SEQ_LENGTH};
            
            if ($pb->is_mapped($flag)) {
                $this_rg_data[2]++;
                $this_rg_data[3] += $pr->{MAPPED_SEQ_LENGTH};
                $this_rg_data[4]++ if $pb->is_sequencing_paired($flag);
                
                # avg quality of mapped bases
                foreach my $qual (VRPipe::Parser::fastq->qual_to_ints($pr->{QUAL})) {
                    $this_rg_data[5]++;
                    
                    if ($this_rg_data[5] == 1) {
                        $this_rg_data[6] = $qual;
                    }
                    else {
                        $this_rg_data[6] += ($qual - $this_rg_data[6]) / $this_rg_data[5];
                    }
                }
                
                # avg insert size and keep track of 's' for later calculation of sd.
                # algorithm based on http://www.johndcook.com/standard_deviation.html
                if ($pb->is_mapped_paired($flag)) {
                    $this_rg_data[7]++;
                    
                    my $isize = $pr->{ISIZE} || 0;
                    
                    if ($pr->{MAPQ} > 0) {
                        if ($isize > 0) { # avoids counting the isize twice for a pair, since one will be negative
                            $this_rg_data[8]++;
                            
                            if ($this_rg_data[8] == 1) {
                                $this_rg_data[9]  = $isize;
                                $this_rg_data[10] = 0;
                            }
                            else {
                                my $old_mean = $this_rg_data[9];
                                $this_rg_data[9] += ($isize - $old_mean) / $this_rg_data[8];
                                $this_rg_data[10] += ($isize - $old_mean) * ($isize - $this_rg_data[9]);
                            }
                            
                            # also, median insert size. Couldn't find an
                            # accurate running algorithm, but just keeping a
                            # histogram is accurate and uses less than 1MB. We
                            # can use the same histogram later to calculate the
                            # MAD.
                            $this_rg_data[11]->{$isize}++;
                        }
                    }
                }
                
                # for later calculation of mismatch %
                my $nm = $pr->{NM};
                if (defined $nm && $nm ne '*') {
                    $this_rg_data[12] += $pr->{SEQ_LENGTH};
                    $this_rg_data[13] += $nm;
                }
            }
            
            if ($pb->is_duplicate($flag)) {
                $this_rg_data[14]++;
                $this_rg_data[15] += $pr->{SEQ_LENGTH};
            }
            
            $readgroup_data{$rg} = \@this_rg_data;
        }
        
        # calculate the means etc.
        my %stats;
        my $math_util = VRPipe::Utils::Math->new();
        foreach my $rg (sort keys %readgroup_data) {
            my @data = @{ $readgroup_data{$rg} };
            
            # calculate/round stats
            my $avg_qual  = defined $data[5] ? sprintf("%0.2f", $data[6])                   : 0;
            my $avg_isize = defined $data[8] ? sprintf("%0.0f", $data[9])                   : 0;
            my $sd_isize  = $avg_isize       ? sprintf("%0.2f", sqrt($data[10] / $data[8])) : 0;
            my $percent_mismatch = defined $data[13] ? sprintf("%0.2f", (100 / $data[12]) * $data[13]) : 0;
            
            my $median_isize = 0;
            my $mad          = 0;
            if (defined $data[11]) {
                $median_isize = $math_util->histogram_median($data[11]);
                
                my %ads;
                while (my ($isize, $freq) = each %{ $data[11] }) {
                    my $ad = abs($median_isize - $isize);
                    $ads{$ad} += $freq;
                }
                
                $mad = $math_util->histogram_median(\%ads);
            }
            
            my %rg_stats;
            $rg_stats{total_reads}                  = $data[0];
            $rg_stats{total_bases}                  = $data[1];
            $rg_stats{mapped_reads}                 = $data[2];
            $rg_stats{mapped_bases}                 = $data[3];
            $rg_stats{mapped_reads_paired_in_seq}   = $data[4];
            $rg_stats{mapped_reads_properly_paired} = $data[7];
            $rg_stats{avg_qual}                     = $avg_qual;
            $rg_stats{avg_isize}                    = $avg_isize;
            $rg_stats{sd_isize}                     = $sd_isize;
            $rg_stats{percent_mismatch}             = $percent_mismatch;
            $rg_stats{median_isize}                 = $median_isize;
            $rg_stats{mad}                          = $mad;
            $rg_stats{duplicate_reads}              = $data[14] || 0;
            $rg_stats{duplicate_bases}              = $data[15] || 0;
            $stats{$rg}                             = \%rg_stats;
        }
        
        return %stats;
    }

=head2 bas
 
 Title   : bas
 Usage   : $obj->bas('in.bam', 'out.bas');
           $obj->bas('in.bam', 'out.bas', release_date => 'YYYYMMDD',
                                          sequence_index => 'sequence.index',
                                          rg_from_pu => 1);
 Function: Generate a 'bas' file. These are bam statistic files that provide
           various stats about each readgroup in a bam.
 Returns : boolean (true on success)
 Args    : input bam file (must have been run via samtools fillmd so that there
           are accurate NM tags for each record), output filename.
           Optionally, to have DCC-style filenames, supply a release_date =>
           YYYYMMDD int describing the date of the release.
           You can also supply sequence_index => sequence.index file (optional
           if your bam headers are good), and an optional boolean
           rg_from_pu which if true means that RG ids are taken to be
           meaningless, with the true read group identifier being in the PU
           field of the RG line(s) in the header. Finally there is an optional
           chr => Str arguement, which will put that chr in the DCC-style
           filename of column 1.

=cut
    
    method bas (ClassName|Object $self: Str|File $in_bam!, Str|File $out_bas!, Int :$release_date? where {/^\d{8}$/}, Str|File :$sequence_index?, Bool :$rg_from_pu?, Str :$chr?) {
        unless (ref($in_bam) && ref($in_bam) eq 'VRPipe::File') {
            $in_bam = VRPipe::File->create(path => file($in_bam));
        }
        unless (ref($out_bas) && ref($out_bas) eq 'VRPipe::File') {
            $out_bas = VRPipe::File->create(path => file($out_bas));
        }
        
        my $bas_fh = $out_bas->openw;
        
        my $expected_lines = 0;
        
        # print header
        print $bas_fh join("\t", 'bam_filename', 'md5', 'study', 'sample', 'platform', 'library', 'readgroup', '#_total_bases', '#_mapped_bases', '#_total_reads', '#_mapped_reads', '#_mapped_reads_paired_in_sequencing', '#_mapped_reads_properly_paired', '%_of_mismatched_bases', 'average_quality_of_mapped_bases', 'mean_insert_size', 'insert_size_sd', 'median_insert_size', 'insert_size_median_absolute_deviation', '#_duplicate_reads', '#_duplicate_bases'), "\n";
        $expected_lines++;
        
        # get the stats for each read group
        my %readgroup_data = $self->bam_statistics($in_bam->path);
        
        # figure out filename for column 1
        my $dcc_filename;
        if ($release_date) {
            my %readgroup_info;
            my $algorithm = 'unknown_algorithm';
            my $bp = VRPipe::Parser->create('bam', { file => $in_bam });
            %readgroup_info = $bp->readgroup_info();
            $algorithm      = $bp->program;
            $bp->close;
            
            my $sample;
            my $platform = 'unknown_platform';
            my $raw_pl;
            my $project;
            my %techs;
            my $example_rg;
            while (my ($rg, $info) = each %readgroup_info) {
                # there should only be one sample, so we just pick the first
                $sample ||= $info->{SM};
                $example_rg ||= $rg_from_pu ? $info->{PU} : $rg;
                
                # might be more than one of these if we're a sample-level bam.
                # We standardise on the DCC nomenclature for the 3 platforms;
                # they should be in this form anyway, so this is just-in-case
                $platform = $info->{PL};
                $raw_pl   = $platform;
                if ($platform =~ /illumina|slx/i) {
                    $platform = 'ILLUMINA';
                }
                elsif ($platform =~ /solid/i) {
                    $platform = 'ABI_SOLID';
                }
                elsif ($platform =~ /454/) {
                    $platform = 'LS454';
                }
                $techs{$platform}++;
                
                $project ||= $info->{DS};
            }
            $sample ||= 'unknown_sample';
            
            # picard merge may have tried to uniqueify the rg, so pluck off .\d
            $example_rg =~ s/\.\d+$//;
            
            unless (defined $example_rg) {
                $self->throw("The bam '" . $in_bam->path . "' had no RG in the header!");
            }
            
            # instead of the srp (project code), we now have population and
            # analysis group in it's place. These things are not stored in the
            # bam header, so we must check the sequence.index file to figure
            # this out. We assume that all the readgroups have the same pop and
            # ag, since they DCC only do same-sample bams
            # sometimes dcc_filename is used just to make bas files, and
            # sometimes we don't care about the dcc_file being correct, so we
            # don't want to require passing in sequence.index file
            my ($pop, $ag) = ('unknown_population', 'unknown_analysisgroup');
            if ($sequence_index) {
                my $sip = VRPipe::Parser->create('sequence_index', { file => $sequence_index });
                $pop = $sip->lane_info($example_rg, 'POPULATION')     || 'unknown_population';
                $ag  = $sip->lane_info($example_rg, 'ANALYSIS_GROUP') || 'unknown_analysisgroup';
                $ag =~ s/\s/_/g;
            }
            
            # if there's more than 1 tech, it doesn't appear in the filename
            if (keys %techs > 1) {
                $platform = '';
            }
            else {
                $platform .= '.';
            }
            
            my $bamname = $in_bam->basename;
            my $chrom   = '';
            if ($chr) {
                if ($chr =~ /^(?:\d+|[XY]|MT)$/) {
                    $chr = 'chrom' . $chr;
                }
                $chrom = "$chr.";
            }
            elsif ($bamname =~ /\.?(chrom(?:\d+|[XY]|MT)|nonchrom|unmapped|mapped)\./) {
                $chrom = "$1.";
            }
            
            # picard merge can fuck with program names, converting them to unique numbers
            if ($algorithm eq 'unknown_algorithm' || $algorithm =~ /^\d+$/ || $algorithm =~ /GATK/) {
                if ($platform =~ /ILLUMINA/) {
                    $algorithm = 'bwa';
                }
                elsif ($platform =~ /454/) {
                    $algorithm = 'ssaha2';
                }
                elsif ($raw_pl eq 'solid') {
                    $algorithm = 'mosaik';
                }
                elsif ($platform =~ /SOLID/) {
                    $algorithm = 'bfast';
                }
            }
            
            $dcc_filename = "$sample.$chrom$platform$algorithm.$pop.$ag.$release_date";
        }
        else {
            $dcc_filename = $in_bam->basename;
        }
        
        # get the md5 of the bam file
        my $md5 = $in_bam->md5;
        unless ($md5) {
            $in_bam->update_md5;
            $md5 = $in_bam->md5;
        }
        
        my $sip = VRPipe::Parser->create('sequence_index', { file => $sequence_index }) if $sequence_index;
        my $pb = VRPipe::Parser->create('bam', { file => $in_bam });
        
        # if necessary, convert from arbitrary RG ids to the lane identifier stored
        # in PU
        my %orig_rgs;
        if ($rg_from_pu) {
            while (my ($rg, $data) = each %readgroup_data) {
                my $new_rg = $pb->readgroup_info($rg, 'PU');
                $orig_rgs{$new_rg} = $rg;
            }
            while (my ($new_rg, $old_rg) = each %orig_rgs) {
                my $data = delete $readgroup_data{$old_rg};
                $readgroup_data{$new_rg} = $data;
            }
        }
        else {
            %orig_rgs = map { $_ => $_ } keys %readgroup_data;
        }
        
        # add in the meta data
        while (my ($rg, $data) = each %readgroup_data) {
            $readgroup_data{$rg}->{dcc_filename} = $dcc_filename;
            $readgroup_data{$rg}->{study}        = $pb->readgroup_info($orig_rgs{$rg}, 'DS') || 'unknown_study';
            $readgroup_data{$rg}->{sample}       = $pb->readgroup_info($orig_rgs{$rg}, 'SM') || 'unknown_sample';
            $readgroup_data{$rg}->{platform}     = $pb->readgroup_info($orig_rgs{$rg}, 'PL') || 'unknown_platform';
            $readgroup_data{$rg}->{library}      = $pb->readgroup_info($orig_rgs{$rg}, 'LB') || 'unknown_library';
            
            # fall back on the sequence.index if we have to
            if ($sip) {
                my %convert = (study    => 'STUDY_ID',
                               sample   => 'SAMPLE_NAME',
                               platform => 'INSTRUMENT_PLATFORM',
                               library  => 'LIBRARY_NAME');
                foreach my $field (keys %convert) {
                    if (!$readgroup_data{$rg}->{$field} || $readgroup_data{$rg}->{$field} eq "unknown_$field" || $readgroup_data{$rg}->{$field} eq '-') {
                        $readgroup_data{$rg}->{$field} = $sip->lane_info($rg, $convert{$field}) || "unknown_$field";
                    }
                }
            }
            
            if (defined $tech_to_platform{ $readgroup_data{$rg}->{platform} }) {
                $readgroup_data{$rg}->{platform} = $tech_to_platform{ $readgroup_data{$rg}->{platform} };
            }
            
            $readgroup_data{$rg}->{md5} = $md5;
        }
        
        # print stats per-readgroup
        foreach my $rg (sort keys %readgroup_data) {
            my %data = %{ $readgroup_data{$rg} };
            print $bas_fh join("\t", $data{dcc_filename}, $data{md5}, $data{study}, $data{sample}, $data{platform}, $data{library}, $rg, $data{total_bases} || 0, $data{mapped_bases} || 0, $data{total_reads} || 0, $data{mapped_reads} || 0, $data{mapped_reads_paired_in_seq} || 0, $data{mapped_reads_properly_paired} || 0, $data{percent_mismatch}, $data{avg_qual}, $data{avg_isize}, $data{sd_isize}, $data{median_isize}, $data{mad}, $data{duplicate_reads} || 0, $data{duplicate_bases} || 0), "\n";
            $expected_lines++;
        }
        $out_bas->close;
        
        my $actual_lines = $out_bas->lines;
        if ($actual_lines == $expected_lines) {
            return 1;
        }
        else {
            $out_bas->unlink;
            $self->warn("Wrote $expected_lines lines to " . $out_bas->path . ", but only read back $actual_lines! Deleted the output.");
            return 0;
        }
    }
}

1;
