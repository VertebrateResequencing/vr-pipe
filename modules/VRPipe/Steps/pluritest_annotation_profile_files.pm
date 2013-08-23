
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

class VRPipe::Steps::pluritest_annotation_profile_files  with VRPipe::StepRole  {
    use File::Compare qw(compare_text);
    
    method options_definition {
        return {
            annot_file_regex   => VRPipe::StepOption->create(description => 'regex used to search for annotation files amongst the expression analysis files obtained from iRODS',     default_value => 'annotation.txt'),
            profile_file_regex => VRPipe::StepOption->create(description => 'regex used to search for sample profile files amongst the expression analysis files obtained from iRODS', default_value => 'Sample_Probe_Profile.txt'),
            header_regex       => VRPipe::StepOption->create(description => 'regex used to identify the header line of the annotation file',                                           default_value => '^TargetID'),
            control_regex      => VRPipe::StepOption->create(description => 'regex used to identify a control sample/lane from the idat file metadata',                                default_value => 'Control'),
        };
    }
    
    method inputs_definition {
        return {
            idat_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'idat file with associated genome studio irods file path and analysis_uuid metadata',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' }
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
                max_files   => 1,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
            ),
            #~ merged_profile_file => VRPipe::StepIODefinition->create(
            #~ type        => 'txt',
            #~ description => 'a file that contains the GenomeStudio profile for the samples',
            #~ max_files   => 1,
            #~ ),
        };
    }
    
    method body_sub {
        return sub {
            my $self               = shift;
            my $options            = $self->options;
            my $annot_file_regex   = $options->{annot_file_regex};
            my $profile_file_regex = $options->{profile_file_regex};
            my $header_regex       = $options->{header_regex};
            my $control_regex      = $options->{control_regex};
            my $req                = $self->new_requirements(memory => 500, time => 1);
            
            my %annotation_files;
            my %profile_files;
            my @annot_profile_map;
            my @profile_lanes;
            my @sample_tag_map;
            
            foreach my $idat_file (@{ $self->inputs->{idat_files} }) {
                my $meta         = $idat_file->metadata;
                my $storage_path = $meta->{storage_path};
                my $control      = $meta->{control};
                my $lib_tag      = $meta->{library_tag};
                my $sample       = $meta->{sample};
                $sample = $sample . "_CTRL" if $control eq $control_regex;
                my $sample_tag = (split /_/, $lib_tag)[0] . '_' . $sample;
                push @sample_tag_map, ($sample_tag, $lib_tag);
                push @profile_lanes, $meta->{lane};
                my @expression_files = <$storage_path/*>;
                my ($annot_file, $pro_file);
                
                foreach my $annotation_path (@expression_files) {
                    $annot_file = $annotation_path if $annotation_path =~ m/$annot_file_regex$/;
                    $pro_file   = $annotation_path if $annotation_path =~ m/$profile_file_regex$/;
                }
                if ($annot_file && !$annotation_files{$annot_file}) {
                    $annotation_files{$annot_file} = 1;
                    push @annot_profile_map, ($annot_file, $pro_file);
                }
                push @{ $profile_files{$pro_file} }, $lib_tag if $pro_file && !$profile_files{$pro_file};
            }
            my @annotation_paths = keys %annotation_files;
            my @profile_paths    = keys %profile_files;
            my @chars            = ("A" .. "Z");
            my $random;
            $random .= $chars[rand @chars] for 1 .. 8;
            my $merged_annotation_file      = $self->output_file(output_key => 'annotation_file', basename => $random . '_annotation.txt', type => 'txt');
            my $merged_annotation_file_path = $merged_annotation_file->path;
            my $mapping_file                = $self->output_file(output_key => 'mapping_file', basename => $random . '_mapping.txt', type => 'txt');
            my $mapping_file_path           = $mapping_file->path;
            my $merged_profile_file         = $self->output_file(temporary => 1, basename => $random . '_merged_profile.txt', type => 'txt');
            my $merged_profile_file_path    = $merged_profile_file->path;
            my $profile_file                = $self->output_file(output_key => 'profile_file', basename => $random . '_profile.txt', type => 'txt');
            my $profile_file_path           = $profile_file->path;
            my $this_cmd                    = "use VRPipe::Steps::pluritest_annotation_profile_files; VRPipe::Steps::pluritest_annotation_profile_files->merge_annotation_files(q[$merged_annotation_file_path], q[$mapping_file_path], q[$profile_file_path], q[$merged_profile_file_path], q[$header_regex], q[$random], profile_paths => [qw(@profile_paths)], annotation_paths => [qw(@annotation_paths)], sample_tags => [qw(@sample_tag_map)], profile_lanes => [qw(@profile_lanes)], annot_profile_map => [qw(@annot_profile_map)]);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$merged_annotation_file, $mapping_file, $profile_file, $merged_profile_file] });
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
    
    method merge_annotation_files (ClassName|Object $self: Str|File $merged_out!, Str|File $mapping_out!, Str|File $profile_out!, Str|File $merged_profile_out!, Str $header_regex!, Str $random!, ArrayRef[Str] :$profile_paths!, ArrayRef[Str] :$annotation_paths!, ArrayRef[Str] :$sample_tags!, ArrayRef[Str] :$profile_lanes!, ArrayRef[Str] :$annot_profile_map!) {
        my %annot_profiles  = @{$annot_profile_map};
        my %sample_tag_hash = @{$sample_tags};
        my @samples         = keys %sample_tag_hash;
        @samples > 0 || $self->throw("No samples and tags provided for mapping file, so can not proceed\n");
        my @expression_order;
        my $new_meta = { merge_tag_id => $random, lanes => join(',', @{$profile_lanes}) };
        my $mapping_file = VRPipe::File->get(path => $mapping_out);
        my $map_fh = $mapping_file->openw;
        while (my ($sample, $tag) = each %sample_tag_hash) {
            push @expression_order, $tag;
            print $map_fh "$tag\t$sample\n";
        }
        $mapping_file->update_stats_from_disc(retries => 3);
        $mapping_file->add_metadata($new_meta);
        $mapping_file->close;
        
        @expression_order = sort @expression_order;
        my @annot_paths = keys %annot_profiles;
        
        system("rm $merged_profile_out");
        my $profile_file = VRPipe::File->get(path => $profile_out);
        my $new_profile = 1;
        
        @annot_paths > 0 || $self->throw("No annotation files provided for the merge.");
        while (@annot_paths > 0) {
            my $annot_path = pop @annot_paths;
            if (!-f $merged_out) { #No file, so create one....
                my $merged_file = VRPipe::File->get(path => $merged_out);
                my $merge_fh = $merged_file->openw;
                my $fh;
                open($fh, $annot_path) || $self->throw("Couldn't open '$annot_path': $!");
                my $body_text = 0;
                while (<$fh>) {
                    $body_text = 1 if $_ =~ /$header_regex/;
                    print $merge_fh $_ if $body_text;
                }
                $merged_file->update_stats_from_disc(retries => 3);
                $merged_file->add_metadata($new_meta);
                $merged_file->close;
            }
            elsif (compare_text($merged_out, $annot_path)) { #Files not identical, so perform intersection
                my %current_fields;
                my @new_fields;
                my $header_line;
                
                my $cfh;
                open($cfh, $merged_out) || $self->throw("Couldn't open '$merged_out': $!");
                while (<$cfh>) {
                    if ($_ =~ /$header_regex/) {
                        $header_line = $_;
                    }
                    else {
                        $current_fields{$_} = 1;
                    }
                }
                close $cfh;
                
                my $nfh;
                open($nfh, $annot_path) || $self->throw("Couldn't open '$annot_path': $!");
                while (<$nfh>) {
                    push @new_fields, $_ unless $_ =~ /$header_regex/;
                }
                close $nfh;
                
                # the intersection of current and new:
                my @intersection = grep($current_fields{$_}, @new_fields);
                my $merged_file = VRPipe::File->get(path => $merged_out);
                my $merge_fh = $merged_file->openw;
                print $merge_fh $header_line;
                for my $annot (@intersection) {
                    print $merge_fh $annot;
                }
                $merged_file->update_stats_from_disc(retries => 3);
                $merged_file->add_metadata($new_meta);
                $merged_file->close;
            }
            
            my $profile_in = $annot_profiles{$annot_path} || $self->throw("Couldn't open '$annot_profiles{$annot_path}': $!");
            if (!-f $merged_profile_out) { #No file, so create one....
                my $merged_profile_file = VRPipe::File->get(path => $merged_profile_out);
                my $pro_fh = $merged_profile_file->openw;
                my $pfh;
                open($pfh, $profile_in) || $self->throw("Couldn't open '$profile_in': $!");
                my $body_text = 0;
                while (<$pfh>) {
                    $body_text = 1 if $_ =~ /$header_regex/;
                    print $pro_fh $_ if $body_text;
                }
            }
            elsif (compare_text($merged_profile_out, $profile_in)) {
                #Files not identical, so perform intersection
                my %current_fields;
                my %new_profiles;
                my @new_fields;
                my @header_line;
                my $header;
                
                my $cfh;
                open($cfh, $merged_profile_out) || $self->throw("Couldn't open '$merged_profile_out': $!");
                while (<$cfh>) {
                    chomp;
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
                
                my $nfh;
                open($nfh, $profile_in) || $self->throw("Couldn't open '$profile_in': $!");
                my $body_text = 0;
                while (<$nfh>) {
                    chomp;
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
                
                my $merged_profile_file = VRPipe::File->get(path => $merged_profile_out);
                my $mp_fh = $merged_profile_file->openw;
                print $mp_fh $header;
                for my $annot (@new_fields) {
                    if (defined $current_fields{$annot} && defined $new_profiles{$annot}) {
                        my $line = $current_fields{$annot} . "\t" . $new_profiles{$annot} . "\n";
                        print $mp_fh $line;
                    }
                }
            }
        }
        
        my $pfh;
        open($pfh, $merged_profile_out) || $self->throw("Couldn't open '$merged_profile_out': $!");
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
        $profile_file->add_metadata($new_meta);
        $profile_file->close;
        return 1;
    }
}

1;
