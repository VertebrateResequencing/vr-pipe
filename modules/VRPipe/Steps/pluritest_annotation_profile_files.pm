
=head1 NAME

VRPipe::Steps::pluritest_annotation_profile_files - a step

=head1 DESCRIPTION

This step produces an annotation file for pluritest analyses that contains an
intersection of  the annotation files provided for the samples/cohorts (i.e.
the annotation file produce contains only those gene expression targets that
are common to all annotation files provided here).

It also produces a mapping file containing a map of sample name to 'analysis
tag', which is a letter corresponding to  the position of the sample data in
the Sample_Probe_Profile file, e.g. 'I' => 'qc1hip5529783', 'K' =>
'qc1hip5529784'.

The third file produced is the profile file containing only those profiles for
the samples involved, which is extracted  from the large Sample_Probe_Profile
file.

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013,2014 Genome Research Limited.

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

class VRPipe::Steps::pluritest_annotation_profile_files  with VRPipe::StepRole  {
    use File::Compare qw(compare_text);
    
    method options_definition {
        return {
            profile_normalisation_method => VRPipe::StepOption->create(
                description   => 'the normalisation method used for generating the sample probe profile file',
                default_value => 'none'
            ),
            header_regex => VRPipe::StepOption->create(
                description   => 'regex used to identify the header line of the annotation file',
                default_value => '^TargetID'
            ),
            mapping_sample_from_metadata => VRPipe::StepOption->create(
                description   => 'metadata key from which the sample name will be taken, for use in the mapping file; separate multiple keys with + symbols - values will be joined with underscores',
                optional      => 1,
                default_value => 'public_name+sample'
            ),
        };
    }
    
    method inputs_definition {
        return {
            idat_files => VRPipe::StepIODefinition->create(
                type            => 'idat',
                description     => 'idat files for samples of the same individual, with associated genome studio irods metadata, where analysis files have already been downloaded',
                max_files       => -1,
                check_existence => 0,
                metadata        => { sample => 'sample name for cell line', irods_analysis_files => 'full irods path to genome studio analysis files', irods_local_storage_dir => 'local base directory where the irods_analysis_files were downloaded', analysis_uuid => 'analysis_uuid', beadchip => 'the chip identifer', beadchip_section => 'section of the beadchip', sample_control => 'boolean for if this sample is a control' }
            ),
        };
    }
    
    method outputs_definition {
        return {
            annotation_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'the annotation file that is used by pluritest - if more than one is provided, the intersection of the annotation files is produced ',
                max_files   => 1
            ),
            mapping_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that maps the samples to the columns in the sample profile file',
                max_files   => 1
            ),
            profile_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that contains the GenomeStudio profile for the samples',
                max_files   => 1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self         = shift;
            my $options      = $self->options;
            my $profile_nm   = $options->{profile_normalisation_method};
            my $header_regex = $options->{header_regex};
            my $mfm          = $options->{mapping_sample_from_metadata};
            my @mfm_keys     = split(/\+/, $mfm);
            my $req          = $self->new_requirements(memory => 500, time => 1);
            
            my %annotation_files;
            my %profile_files;
            my @annot_profile_map;
            my @sample_tag_map;
            
            my $merged_meta = $self->combined_metadata($self->inputs->{idat_files});
            
            my $merged_annotation_file      = $self->output_file(output_key => 'annotation_file', basename => 'annotation.txt', type => 'txt', metadata => $merged_meta);
            my $merged_annotation_file_path = $merged_annotation_file->path->stringify;
            my $mapping_file                = $self->output_file(output_key => 'mapping_file', basename => 'mapping.txt', type => 'txt', metadata => $merged_meta);
            my $mapping_file_path           = $mapping_file->path->stringify;
            my $profile_file                = $self->output_file(output_key => 'profile_file', basename => 'profile.txt', type => 'txt', metadata => $merged_meta);
            my $profile_file_path           = $profile_file->path->stringify;
            
            foreach my $idat_file (@{ $self->inputs->{idat_files} }) {
                my $idat_path = $idat_file->path->stringify;
                $self->relate_input_to_output($idat_path, 'individual_annotation_merge', $merged_annotation_file_path);
                $self->relate_input_to_output($idat_path, 'individual_mapping_merge',    $mapping_file_path);
                $self->relate_input_to_output($idat_path, 'individual_profile_merge',    $profile_file_path);
                
                my $meta           = $idat_file->metadata;
                my $analysis_files = $meta->{irods_analysis_files};
                my $local_dir      = $meta->{irods_local_storage_dir};
                my $lib_tag        = $meta->{beadchip} . '_' . $meta->{beadchip_section};
                my $sample         = $meta->{sample};
                my $control        = $meta->{sample_control};
                my $sample_tag     = join('_', map { $meta->{$_} } @mfm_keys);
                $sample_tag .= '_CTRL' if $control;
                push @sample_tag_map, ($sample_tag, $lib_tag);
                
                my ($annot_file, $pro_file);
                foreach my $path (@$analysis_files) {
                    $path = file($local_dir, $path)->stringify;
                    my $file = VRPipe::File->get(path => $path);
                    my $meta = $file->metadata;
                    my $type = $meta->{summary_type} || next;
                    if ($type eq 'annotation') {
                        $annot_file = $path;
                    }
                    elsif ($type eq 'probe') {
                        my $group = $meta->{summary_group} || next;
                        next unless $group eq 'sample';
                        my $this_nm = $meta->{normalisation_method} || next;
                        $pro_file = $path if $this_nm eq $profile_nm;
                    }
                }
                
                unless ($annot_file && $pro_file) {
                    $self->throw("idat file $idat_path did not have both an annotation file and suitable sample probe profile file associated with it");
                }
                
                if ($annot_file && !$annotation_files{$annot_file}) {
                    $annotation_files{$annot_file} = 1;
                    push @annot_profile_map, ($annot_file, $pro_file);
                }
                push @{ $profile_files{$pro_file} }, $lib_tag if $pro_file && !$profile_files{$pro_file};
            }
            
            my $merged_profile_file = $self->output_file(temporary => 1, basename => 'merged_profile.txt', type => 'txt');
            my $merged_profile_file_path = $merged_profile_file->path;
            
            my @annotation_paths = keys %annotation_files;
            my @profile_paths    = keys %profile_files;
            my $this_cmd         = "use VRPipe::Steps::pluritest_annotation_profile_files; VRPipe::Steps::pluritest_annotation_profile_files->merge_annotation_files(q[$merged_annotation_file_path], q[$mapping_file_path], q[$profile_file_path], q[$merged_profile_file_path], q[$header_regex], profile_paths => [qw(@profile_paths)], annotation_paths => [qw(@annotation_paths)], sample_tags => [qw(@sample_tag_map)], annot_profile_map => [qw(@annot_profile_map)]);";
            $self->dispatch_vrpipecode($this_cmd, $req);
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method description {
        return "This step produces a merged annotation file as well as mapping and profile files for the selected samples for pluritest analyses.";
    }
    
    method merge_annotation_files (ClassName|Object $self: Str|File $merged_out!, Str|File $mapping_out!, Str|File $profile_out!, Str|File $merged_profile_out!, Str $header_regex!, ArrayRef[Str] :$profile_paths!, ArrayRef[Str] :$annotation_paths!, ArrayRef[Str] :$sample_tags!, ArrayRef[Str] :$annot_profile_map!) {
        my %annot_profiles  = @{$annot_profile_map};
        my %sample_tag_hash = @{$sample_tags};
        my @samples         = keys %sample_tag_hash;
        @samples > 0 || $self->throw("No samples and tags provided for mapping file, so can not proceed\n");
        my @expression_order;
        my $mapping_file = VRPipe::File->get(path => $mapping_out);
        my $map_fh = $mapping_file->openw;
        while (my ($sample, $tag) = each %sample_tag_hash) {
            push @expression_order, $tag;
            print $map_fh "$tag\t$sample\n";
        }
        $mapping_file->update_stats_from_disc(retries => 3);
        $mapping_file->close;
        
        @expression_order = sort @expression_order;
        my @all_annot_paths = keys %annot_profiles;
        
        @all_annot_paths > 0 || $self->throw("No annotation files provided for the merge.");
        
        my @req_cols = ("TargetID", "ProbeID", "REFSEQ_ID", "SYMBOL", "PROBE_ID", "PROBE_SEQUENCE", "CHROMOSOME", "PROBE_COORDINATES");
        
        my %miss_cols = ();
        foreach my $annot_path (@all_annot_paths) {
            open(my $fh, $annot_path) || die "Could not open $annot_path\n";
            while (<$fh>) {
                next unless $_ =~ /$header_regex/;
                for my $col (0 .. $#req_cols) {
                    my $pattern = $req_cols[$col];
                    $miss_cols{$col} = $annot_path unless $_ =~ qr/$pattern/;
                }
                last;
            }
            close($fh);
        }
        
        my $merged_profile_file = VRPipe::File->get(path => $merged_profile_out);
        $merged_profile_file->remove;
        my $profile_file = VRPipe::File->get(path => $profile_out);
        my $merged_file  = VRPipe::File->get(path => $merged_out);
        my $new_profile  = 1;
        
        my $header_line = join("\t", @req_cols);
        while (@all_annot_paths > 0) {
            my $annot_path = pop @all_annot_paths;
            my $annot_file = VRPipe::File->get(path => $annot_path);
            
            if (!$merged_file->e) { #No file, so create one....
                my $merge_fh  = $merged_file->openw;
                my $fh        = $annot_file->openr;
                my $body_text = 0;
                while (my $line = <$fh>) {
                    if ($line =~ /$header_regex/) {
                        print $merge_fh $header_line . "\n";
                        $body_text = 1;
                    }
                    elsif ($body_text) {
                        $line =~ s/\t\r$//;
                        my @line_arr = split(/\t/, $line);
                        my @new_line = ();
                        my $j        = 0;
                        for my $i (0 .. $#req_cols) {
                            if (exists $miss_cols{$i}) {
                                push(@new_line, "NaN");
                                $j++ if ($miss_cols{$i} eq $annot_path);
                            }
                            else {
                                push(@new_line, $line_arr[$i - $j]);
                            }
                        }
                        print $merge_fh join("\t", @new_line);
                    }
                }
                $merged_file->update_stats_from_disc(retries => 3);
                $merged_file->close;
            }
            elsif (compare_text($merged_out, $annot_path)) { #Files not identical, so perform intersection
                my %current_fields;
                my @new_fields = ();
                
                my $cfh = $merged_file->openr;
                while (<$cfh>) {
                    $current_fields{$_} = 1;
                }
                close $cfh;
                
                my $nfh       = $annot_file->openr;
                my $body_text = 0;
                while (my $line = <$nfh>) {
                    if ($line =~ /$header_regex/) {
                        $body_text = 1;
                    }
                    elsif ($body_text) {
                        $line =~ s/\t\r$//;
                        my @line_arr = split(/\t/, $line);
                        my @new_line = ();
                        my $j        = 0;
                        for my $i (0 .. $#req_cols) {
                            if (exists $miss_cols{$i}) {
                                push(@new_line, "NaN");
                                $j++ if ($miss_cols{$i} eq $annot_path);
                            }
                            else {
                                push(@new_line, $line_arr[$i - $j]);
                            }
                        }
                        my $new_field = join("\t", @new_line);
                        push(@new_fields, $new_field);
                    }
                }
                $annot_file->close;
                
                # the intersection of current and new:
                my @intersection = grep($current_fields{$_}, @new_fields);
                my $merge_fh = $merged_file->openw;
                print $merge_fh $header_line . "\n";
                for my $annot (@intersection) {
                    print $merge_fh $annot;
                }
                $merged_file->update_stats_from_disc(retries => 3);
                $merged_file->close;
            }
            
            my $profile_in = $annot_profiles{$annot_path} || $self->throw("Couldn't open '$annot_profiles{$annot_path}': $!");
            my $profile_in_file = VRPipe::File->get(path => $profile_in);
            if (!$merged_profile_file->e) { #No file, so create one....
                my $pro_fh    = $merged_profile_file->openw;
                my $pfh       = $profile_in_file->openr;
                my $body_text = 0;
                while (<$pfh>) {
                    chomp;
                    $_ =~ s/\t?\r?$//;
                    $body_text = 1 if $_ =~ /$header_regex/;
                    print $pro_fh $_ . "\n" if $body_text;
                }
            }
            elsif (compare_text($merged_profile_out, $profile_in)) {
                #Files not identical, so perform intersection
                my %current_fields;
                my %new_profiles;
                my @new_fields;
                my @header_line;
                my $header;
                
                my $cfh = $merged_profile_file->openr;
                while (<$cfh>) {
                    chomp;
                    $_ =~ s/\t?\r?$//;
                    my @fields = split /\s+/, $_;
                    if ($_ =~ /$header_regex/) {
                        @header_line = @fields;
                    }
                    else {
                        my $key = "$fields[0]_$fields[1]";
                        $current_fields{$key} = join("\t", @fields);
                    }
                }
                close $cfh;
                
                my $nfh       = $profile_in_file->openr;
                my $body_text = 0;
                while (<$nfh>) {
                    chomp;
                    $_ =~ s/\t?\r?$//;
                    my @fields = split /\s+/, $_;
                    if ($_ =~ /$header_regex/) {
                        $header = join("\t", (@header_line, @fields[2 .. $#fields]));
                        $body_text = 1;
                    }
                    else {
                        if ($body_text) {
                            my $key = "$fields[0]_$fields[1]";
                            push @new_fields, $key;
                            $new_profiles{$key} = join("\t", @fields[2 .. $#fields]);
                        }
                    }
                }
                close $nfh;
                
                my $mp_fh = $merged_profile_file->openw;
                print $mp_fh $header . "\n";
                for my $annot (@new_fields) {
                    if (defined $current_fields{$annot} && defined $new_profiles{$annot}) {
                        my $line = $current_fields{$annot} . "\t" . $new_profiles{$annot} . "\n";
                        print $mp_fh $line;
                    }
                }
            }
        }
        
        my $pfh       = $merged_profile_file->openr;
        my $pro_fh    = $profile_file->openw;
        my $body_text = 0;
        my @profile_order;
        my @head_line;
        while (<$pfh>) {
            chomp;
            my @line = split /\s+/, $_;
            if ($_ =~ /$header_regex/) {
                $body_text = 1;
                @head_line = @line;
                for my $analysis_tag (@expression_order) {
                    for (my $i = 2; $i <= $#head_line; $i++) {
                        if ($head_line[$i] =~ m/$analysis_tag/) {
                            push @profile_order, $i;
                        }
                    }
                }
            }
            if ($body_text) {
                print $pro_fh "$line[0]\t$line[1]";
                for my $key (@profile_order) {
                    print $pro_fh "\t$line[$key]";
                }
                print $pro_fh "\n";
            }
        }
        close $pfh;
        $profile_file->update_stats_from_disc(retries => 3);
        $profile_file->close;
        return 1;
    }
}

1;
