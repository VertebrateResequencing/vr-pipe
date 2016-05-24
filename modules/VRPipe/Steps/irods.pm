
=head1 NAME

VRPipe::Steps::irods - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012,2014-2016 Genome Research Limited.

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

class VRPipe::Steps::irods with VRPipe::StepRole {
    use VRPipe::Schema;
    
    our $schema;
    our $schema_vrtrack;
    
    has 'irods_exes' => (
        is      => 'ro',
        isa     => 'HashRef',
        lazy    => 1,
        builder => '_build_irods_exes'
    );
    
    method _build_irods_exes {
        return {};
    }
    
    method handle_exes (HashRef $options) {
        my $irods_exes = $self->irods_exes;
        foreach my $exe (keys %{$irods_exes}) {
            my $method = $exe . '_exe';
            next unless defined $options->{$method};
            $irods_exes->{$exe} = $options->{$method};
        }
    }
    
    method options_definition {
        my %opts;
        
        my $irods_exes = $self->irods_exes;
        foreach my $exe (keys %{$irods_exes}) {
            $opts{ $exe . '_exe' } = VRPipe::StepOption->create(description => "path to your irods '$exe' executable", optional => 1, default_value => $irods_exes->{$exe} || $exe);
        }
        
        return \%opts;
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub { return 1; };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generic step for steps using irods";
    }
    
    method max_simultaneous {
        return 15;
    }
    
    method get_file_by_basename (ClassName|Object $self: Str :$basename!, Str|File :$dest!, Str :$zone = 'seq', Str|File :$iget!, Str|File :$iquest!, Str|File :$ichksum!, Str|File :$samtools_for_cram_to_bam?, Str :$iget_args?, Str :$ichksum_args?) {
        my $dest_file = VRPipe::File->get(path => $dest);
        $dest_file->disconnect;
        
        my @cmd_output = $self->open_irods_command(qq[$iquest -z $zone "SELECT COLL_NAME, DATA_NAME WHERE DATA_NAME = '$basename'"]);
        my ($path, $filename);
        foreach (@cmd_output) {
            # output looks like:
            # Zone is seq
            # COLL_NAME = /seq/5150
            # DATA_NAME = 5150_1#1.bam
            # -------------------------
            #
            # or an error is thrown, and only:
            # Zone is seq
            # is output
            
            if (/^COLL_NAME = (.+)$/) {
                $path = $1;
            }
            if (/^DATA_NAME = (.+)$/) {
                $filename = $1;
            }
        }
        
        if ($path && $filename) {
            unless ($filename eq $basename) {
                $self->throw("$filename should be the same as $basename");
            }
            my $irods_file = join('/', ($path, $filename));
            
            $self->get_file(source => $irods_file, dest => $dest_file->path, iget => $iget, ichksum => $ichksum, $samtools_for_cram_to_bam ? (samtools_for_cram_to_bam => $samtools_for_cram_to_bam) : (), $iget_args ? (iget_args => $iget_args) : (), $ichksum_args ? (ichksum_args => $ichksum_args) : ());
        }
        else {
            $self->throw("A file with basename $basename could not be found in iRods zone $zone");
        }
    }
    
    method get_file_md5 (ClassName|Object $self: Str :$file!, Str|File :$ichksum!, Str :$ichksum_args!) {
        $file = $self->_path_escape($file);
        my ($chksum) = $self->open_irods_command("$ichksum $ichksum_args $file");
        my ($md5)    = $chksum =~ m/\b([0-9a-f]{32})\b/i;                        # 32 char MD5 hex string
        unless ($md5) {
            $self->throw("Could not determine md5 of $file in IRODS ('$ichksum $ichksum_args $file' returned '$chksum'; aborted");
        }
        return $md5;
    }
    
    method get_file (ClassName|Object $self: Str :$source!, Str|File :$dest!, Str|File :$iget!, Str|File :$ichksum!, Bool :$add_metadata?, Str|File :$imeta?, Str|File :$samtools_for_cram_to_bam?, Str :$iget_args = '-K -f', Str :$ichksum_args = '') {
        my $dest_file = VRPipe::File->get(path => $dest);
        $dest_file->disconnect;
        
        # before we go fetch a file, check the md5 matches what we're expecting
        my $irodschksum = $self->get_file_md5(file => $source, ichksum => $ichksum, ichksum_args => $ichksum_args);
        my $expected_md5 = $dest_file->metadata->{expected_md5} || $irodschksum;
        $expected_md5 = $irodschksum if ref($expected_md5);                      # if we have multiple expected_md5 in an array ref, it can't really be the expected md5 of this file
        unless ($irodschksum eq $expected_md5) {
            $dest_file->unlink;
            $self->throw("expected md5 checksum in metadata did not match md5 of $source in IRODS; aborted");
        }
        
        # before we start to fetch the file, populate its metadata, so that if
        # another setup imports the file and sees the existing file it can
        # create a symlink with the correct metadata; for compatability with
        # most vrpipe steps, cram/bam files need a reads metadata that is the
        # 0x900 sequences count. We hope that the total_reads metadata in irods
        # for the file is this count, and force set it here. If irods didn't
        # have this, or it wasn't the 0x900 number, we still have hope of
        # getting the correct value below, when we check the bam stats file
        my $meta = $self->get_file_metadata($source, $imeta ? (imeta => $imeta) : ());
        if (defined $meta->{total_reads}) {
            if (!$add_metadata) {
                $dest_file->add_metadata({ reads => $meta->{total_reads} }, replace_data => 1);
            }
            else {
                $meta->{reads} = $meta->{total_reads};
                $dest_file->add_metadata($meta, replace_data => 1);
            }
        }
        $schema ||= VRPipe::Schema->create('VRPipe');
        my $source_graph_node = $schema->get('File', { path => $source, protocol => 'irods:' }) || $schema->get('File', { path => $source });
        if ($source_graph_node && $dest =~ /\.(?:bam|cram)$/) {
            # if there's a bamstats associated with the source file, add the
            # metadata that many bam/cram-related steps rely on to the dest
            # file
            $schema_vrtrack ||= VRPipe::Schema->create('VRTrack');
            my $qc_meta = $schema_vrtrack->vrtrack_metadata(node => $source_graph_node);
            
            if ($qc_meta) {
                my $graph_meta = {};
                if (defined $qc_meta->{'manual_qc'}) {
                    $graph_meta->{npg_qc_status} = $qc_meta->{'manual_qc'} ? 1 : 0;
                }
                $graph_meta->{lane}     = $qc_meta->{'vrtrack_lane_unique'};
                $graph_meta->{library}  = $qc_meta->{'vrtrack_library_id'};
                $graph_meta->{sample}   = $qc_meta->{'vrtrack_sample_name'};
                $graph_meta->{study}    = $qc_meta->{'vrtrack_study_id'};
                $graph_meta->{study_id} = $qc_meta->{'vrtrack_study_name'};
                
                if (defined $qc_meta->{'vrtrack_bam_stats_total length'}) {
                    $graph_meta->{bases}           = $qc_meta->{'vrtrack_bam_stats_total length'};
                    $graph_meta->{paired}          = $qc_meta->{'vrtrack_bam_stats_reads properly paired'} ? 1 : 0;
                    $graph_meta->{avg_read_length} = $qc_meta->{'vrtrack_bam_stats_average length'};
                    
                    # again, we require the 0x900 stats, but the user may be
                    # tracking a 0xB00 bam stats file, or something else;
                    # only override reads metadata if it's definitely an
                    # 0x900 count
                    if (defined $qc_meta->{'vrtrack_bam_stats_options'} && $qc_meta->{'vrtrack_bam_stats_options'} =~ /0x900\b/) {
                        #*** 'raw total sequences' is total records, 'sequences' is the 0x900 count, but will this be true in samtools 1.3+?
                        $graph_meta->{reads} = $qc_meta->{'vrtrack_bam_stats_sequences'} || $qc_meta->{'vrtrack_bam_stats_raw total sequences'};
                    }
                }
                
                while (my ($key, $val) = each %$graph_meta) {
                    if (defined $val) {
                        $meta->{$key} = $val;
                    }
                }
                
                $dest_file->add_metadata($meta, replace_data => 1);
            }
        }
        
        my $failed;
        my $converted_cram_to_bam = 0;
        my $source_path           = $self->_path_escape($source);
        my $dest_path             = $self->_path_escape("$dest");
        if ($samtools_for_cram_to_bam && $source =~ /\.cram$/ && $dest =~ /\.bam$/) {
            $samtools_for_cram_to_bam = $self->_path_escape("$samtools_for_cram_to_bam");
            $iget_args =~ s/-?[Kf]//g; #*** doesn't support options that take alphanumeric args... hopefully this doesn't come up...
            $iget_args =~ s/^\s+//;
            $iget_args =~ s/\s+$//;
            $failed                = $self->run_irods_command("$iget $iget_args $source_path - | $samtools_for_cram_to_bam view -b - > $dest_path");
            $converted_cram_to_bam = 1;
        }
        else {
            # -K: checksum
            # -Q: use UDP rather than TCP
            # -f: force overwrite
            $failed = $self->run_irods_command("$iget $iget_args $source_path $dest_path");
        }
        
        $dest_file->update_stats_from_disc;
        if ($failed) {
            $dest_file->unlink;
            $self->throw("iget failed for $source -> $dest");
        }
        
        # correct permissions, just in case
        chmod 0664, $dest;
        
        if ($converted_cram_to_bam) {
            # we can't confirm based on md5, so check based on num records
            # (this number may not exist outside Sanger, and within Sanger it
            #  does not include secondary or supplementary reads, so we do no
            #  check if not present and need greater than or equal the
            #  number of reads expected)
            my $actual   = $dest_file->num_records;
            my $expected = $meta->{total_reads};
            if ($expected && ($expected > $actual)) {
                $dest_file->unlink;
                $self->throw("we converted $source -> $dest, but the bam had $actual reads instead of $expected; deleted");
            }
        }
        else {
            # double-check the md5 (iget -K doesn't always work?)
            my $ok = $dest_file->verify_md5(file($dest), $expected_md5);
            unless ($ok) {
                $dest_file->unlink;
                $self->throw("we got $source -> $dest, but the md5 checksum did not match; deleted");
            }
        }
        
        if ($add_metadata) {
            # correct the md5 if we converted the file
            if ($converted_cram_to_bam && defined $meta->{md5}) {
                $dest_file->update_md5();
                $meta->{md5} = $dest_file->md5();
            }
            
            $dest_file->add_metadata($meta, replace_data => 1);
        }
        
        # relate source file to dest file in the graph database
        if ($source_graph_node) {
            $self->relate_input_to_output($source_graph_node, 'imported', $dest_file->path->stringify);
        }
    }
    
    method get_file_metadata (ClassName|Object $self: Str $path!, Str|File :$imeta = 'imeta') {
        $path = $self->_path_escape($path);
        my @cmd_output = $self->open_irods_command("$imeta ls -d $path");
        my $meta       = {};
        my $attribute;
        foreach (@cmd_output) {
            if (/^attribute:\s+(\S+)/) {
                $attribute = $1;
                undef $attribute if $attribute =~ /^dcterms:/;
                undef $attribute if ($attribute && $attribute =~ /_history$/);
            }
            elsif ($attribute && /^value:\s+(.+)$/) {
                if (exists $meta->{$attribute}) {
                    my $previous = $meta->{$attribute};
                    if (ref($previous)) {
                        $meta->{$attribute} = [@$previous, $1];
                    }
                    else {
                        $meta->{$attribute} = [$previous, $1];
                    }
                }
                else {
                    $meta->{$attribute} = $1;
                }
            }
        }
        
        return $meta;
    }
    
    method search_by_metadata (ClassName|Object $self: HashRef :$metadata!, Str|File :$imeta!, Str :$zone = 'seq') {
        my @meta;
        while (my ($key, $val) = each %$metadata) {
            push(@meta, "$key = '$val'");
        }
        my $meta = join(' and ', @meta);
        $meta || $self->throw("No metadata supplied");
        
        my @cmd_output = $self->open_irods_command("$imeta -z $zone qu -d $meta");
        my (@results, $dir);
        foreach (@cmd_output) {
            # output looks like:
            # collection: /seq/fluidigm/multiplexes
            # dataObj: cgp_fluidigm_snp_info_1000Genomes.tsv
            # ----
            # collection: /seq/fluidigm/multiplexes
            # dataObj: ddd_fluidigm_snp_info_1000Genomes.tsv
            #
            # or nothing is found and:
            # No rows found
            
            if (/^collection: (.+)$/) {
                $dir = $1;
            }
            elsif (/^dataObj: (.+)$/) {
                push(@results, "$dir/$1");
            }
        }
        
        return @results;
    }
    
    method open_irods_command (ClassName|Object $self: Str $cmd) {
        my (@output, $error);
        foreach my $i (1 .. 3) {
            my $ok = open(my $irods, "$cmd |");
            unless ($ok) {
                $error = "Failed to open pipe to [$cmd]";
                warn $error, "\n";
                sleep($i);
                next;
            }
            
            while (<$irods>) {
                chomp;
                push(@output, $_);
            }
            
            $ok = close($irods);
            last if $ok;
            
            $error = "Failed to close pipe to [$cmd]";
            warn $error, "\n";
            sleep($i);
            @output = ();
        }
        
        $self->throw($error) if $error;
        return @output;
    }
    
    method run_irods_command (ClassName|Object $self: Str $cmd) {
        my $failed;
        foreach my $i (1 .. 3) {
            $failed = system($cmd);
            if ($failed) {
                sleep($i);
                next;
            }
            last;
        }
        return $failed;
    }
    
    method _path_escape (ClassName|Object $self: Str $path) {
        return quotemeta($path);
        #*** strictly speaking quotemeta isn't meant for this, and
        # String::ShellQuote may be more appropriate, and avoiding the shell
        # in *_irods_commands() would be best, but this should do for our
        # purposes...
    }
}

1;
