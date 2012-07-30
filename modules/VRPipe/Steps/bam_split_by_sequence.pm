
=head1 NAME

VRPipe::Steps::bam_split_by_sequence - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::bam_split_by_sequence with VRPipe::StepRole {
    use VRPipe::Parser;
    use VRPipe::Parser::bam;
    
    method options_definition {
        return { split_bam_ignore        => VRPipe::StepOption->create(description => '\'regex\' to ignore sequences with ids that match the regex, eg. ignore => \'^N[TC]_\d+\' to avoid making splits for contigs',                                                                                                              optional => 1),
                 split_bam_merge         => VRPipe::StepOption->create(description => '{ \'regex\' => \'prefix\' } to put sequences with ids that match the regex into a single file with the corresponding prefix. split_bam_ignore takes precendence.',                                                                          optional => 1, default_value => '{}'),
                 split_bam_include_mate  => VRPipe::StepOption->create(description => 'boolean; if true mates mapping to the split sequence will be included',                                                                                                                                                                     optional => 1, default_value => 0),
                 split_bam_non_chrom     => VRPipe::StepOption->create(description => 'boolean; split_bam_ignore takes precendence: does a merge for sequences that don\'t look chromosomal into a \'nonchrom\' prefixed file; does not make individual bams for each non-chromosomsal seq if split_bam_ignore has not been set',  optional => 1, default_value => 1),
                 split_bam_make_unmapped => VRPipe::StepOption->create(description => 'boolean; a file prefixed with \'unmapped\' is made containing all the unmapped reads. NB: only unmapped pairs or singletons appear in this file; the other split bam files will contain unmapped reads where that read\'s mate was mapped', optional => 1, default_value => 0),
                 split_bam_all_unmapped  => VRPipe::StepOption->create(description => 'boolean; when true, the described behaviour of make_unmapped above changes so that the unmapped file contains all unmapped reads, potentially duplicating reads in different split files',                                                  optional => 1, default_value => 0),
                 split_bam_only          => VRPipe::StepOption->create(description => '\'regex\' to only makes splits for sequences that match the regex. split_bam_non_chr is disabled',                                                                                                                                          optional => 1) };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to be split') };
    }
    
    method body_sub {
        return sub {
            use VRPipe::Steps::bam_split_by_sequence;
            
            my $self    = shift;
            my $options = $self->options;
            $self->throw('split_bam_merge is not properly formed') unless ($options->{split_bam_merge} =~ m/^\s*{.*}\s*$/);
            my $merge = eval "$$options{split_bam_merge}";
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $pb = VRPipe::Parser->create('bam', { file => $bam });
                my %all_sequences = $pb->sequence_info();
                
                my $split_dir       = $self->output_root->absolute->stringify;
                my $bam_path        = $bam->path->stringify;
                my $args            = "q[$bam_path], split_dir => q[$split_dir], ignore => q[$$options{split_bam_ignore}], non_chrom => $$options{split_bam_non_chrom}, merge => $$options{split_bam_merge}, make_unmapped => $$options{split_bam_make_unmapped}, all_unmapped => $$options{split_bam_all_unmapped}, only => q[$$options{split_bam_only}], include_mate => $$options{split_bam_include_mate}";
                my $split_bam_files = VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($bam_path, split_dir => $split_dir, ignore => $$options{split_bam_ignore}, non_chrom => $$options{split_bam_non_chrom}, merge => $merge, make_unmapped => $$options{split_bam_make_unmapped}, all_unmapped => $$options{split_bam_all_unmapped}, only => $$options{split_bam_only}, include_mate => $$options{split_bam_include_mate}, pretend => 1);
                my @output_files;
                while (my ($split, $split_bam_path) = each %$split_bam_files) {
                    my $split_bam = Path::Class::File->new($split_bam_path);
                    my $meta      = $bam->metadata;
                    delete $meta->{reads} if (exists $meta->{reads});
                    delete $meta->{bases} if (exists $meta->{bases});
                    push @output_files,
                      $self->output_file(output_key => 'split_bam_files',
                                         basename   => $split_bam->basename,
                                         type       => 'bam',
                                         metadata   => { %$meta, split_source_bam => $bam_path, split_sequence => $split });
                }
                my $this_cmd = "use VRPipe::Steps::bam_split_by_sequence; VRPipe::Steps::bam_split_by_sequence->split_bam_by_sequence($args);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [@output_files] });
            }
        };
    }
    
    method outputs_definition {
        return { split_bam_files => VRPipe::StepIODefinition->create(type        => 'bam',
                                                                     max_files   => -1,
                                                                     description => 'split bams for each sequence the reads were aligned to',
                                                                     metadata    => {
                                                                                   reads            => 'number of reads in the split bam',
                                                                                   split_source_bam => 'The original bam from which this this split sequence was derived',
                                                                                   split_sequence   => 'The sequence prefix label describing the sequence split out into this bam' }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Splits a BAM file into multiple BAM files, one for each sequence the reads were aligned to";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method split_bam_by_sequence (ClassName|Object $self: Str $bam_file!, Str :$split_dir!, Str :$ignore?, Bool :$non_chrom = 1, HashRef[Str] :$merge = {}, Bool :$make_unmapped = 0, Bool :$all_unmapped = 0, Str :$only?, Bool :$pretend = 0, Bool :$include_mate = 0) {
        unless (ref($bam_file) && ref($bam_file) eq 'VRPipe::File') {
            $bam_file = VRPipe::File->get(path => file($bam_file));
        }
        
        if ($only) {
            $non_chrom = 0;
        }
        
        if ($non_chrom) {
            $merge->{'^(?:N[TC]_\d+|GL\d+|hs\S+)'} = 'nonchrom';
            unless ($ignore) {
                $ignore = '^(?:N[TC]_\d+|GL\d+|hs\S+)';
            }
        }
        
        # find out what sequences there are
        my $pb = VRPipe::Parser->create('bam', { file => $bam_file });
        my %all_sequences = $pb->sequence_info();
        
        # work out the output filenames for each sequence
        my %seq_to_bam;
        my %split_bams;
        foreach my $seq (keys %all_sequences) {
            my @prefixes;
            if ($merge) {
                while (my ($regex, $this_prefix) = each %{$merge}) {
                    if ($seq =~ /$regex/) {
                        push(@prefixes, $this_prefix);
                    }
                }
            }
            unless (scalar(@prefixes) || ($ignore && $seq =~ /$ignore/)) {
                @prefixes = $seq =~ /chr/ ? ($seq) : ('chrom' . $seq);
            }
            
            if ($only) {
                next unless ($seq =~ /$only/);
            }
            
            foreach my $prefix (@prefixes) {
                my $split_bam = VRPipe::File->create(path => file($split_dir, $prefix . '.' . $bam_file->basename));
                $split_bams{$prefix} = $split_bam->path->stringify;
                push(@{ $seq_to_bam{$seq} }, $split_bam->path->stringify);
            }
        }
        
        my @get_fields = ('RNAME');
        my ($unmapped_bam_path, $skip_mate_mapped);
        if ($make_unmapped) {
            my $unmapped_bam = VRPipe::File->create(path => file($split_dir, 'unmapped.' . $bam_file->basename));
            $unmapped_bam_path = $unmapped_bam->path->stringify;
            $split_bams{unmapped} = $unmapped_bam_path;
            push(@get_fields, 'FLAG');
            $skip_mate_mapped = $all_unmapped ? 0 : 1;
        }
        if ($include_mate) {
            push(@get_fields, 'MRNM');
        }
        
        if ($pretend) {
            return \%split_bams;
        }
        
        $bam_file->disconnect;
        
        # stream through the bam, outputting all our desired files
        my $parsed_record = $pb->parsed_record();
        $pb->get_fields(@get_fields);
        my %counts;
        while ($pb->next_record) {
            my $seq_name = $parsed_record->{RNAME};
            my %out_bams;
            if (exists $seq_to_bam{"$seq_name"}) {
                foreach my $out_bam (@{ $seq_to_bam{"$seq_name"} }) {
                    $out_bams{$out_bam} = 1;
                }
            }
            if ($include_mate) {
                my $mate_seq_name = $parsed_record->{MRNM};
                if (exists $seq_to_bam{"$mate_seq_name"}) {
                    foreach my $out_bam (@{ $seq_to_bam{"$mate_seq_name"} }) {
                        $out_bams{$out_bam} = 1;
                    }
                }
            }
            
            if (keys %out_bams) {
                foreach my $out_bam (keys %out_bams) {
                    $pb->write_result($out_bam . '.unchecked.bam');
                    $counts{$out_bam}++;
                }
            }
            
            if ($unmapped_bam_path) {
                my $flag = $parsed_record->{FLAG};
                unless ($pb->is_mapped($flag)) {
                    if ($skip_mate_mapped && $pb->is_sequencing_paired($flag)) {
                        next if $pb->is_mate_mapped($flag);
                    }
                    $pb->write_result($unmapped_bam_path . '.unchecked.bam');
                    $counts{$unmapped_bam_path}++;
                }
            }
        }
        $pb->close;
        
        # check all the bams we created
        my @out_files;
        foreach my $split_path (values %split_bams) {
            my $unchecked_path = $split_path . '.unchecked.bam';
            my $split_bam      = VRPipe::File->get(path => file($split_path));
            my $unchecked      = VRPipe::File->create(path => file($unchecked_path));
            
            $unchecked->update_stats_from_disc(retries => 3);
            
            # the input bam might not have had reads mapped to every sequence in the
            # header, so we might not have created all the bams expected. Create
            # a header-only bam in that case:
            unless ($unchecked->s) {
                $pb = VRPipe::Parser->create('bam', { file => $bam_file });
                $pb->next_record;
                $pb->_create_no_record_output_bam($unchecked_path);
                $pb->close;
            }
            
            $unchecked->disconnect;
            
            my $actual_reads = $unchecked->num_records;
            my $expected_reads = $counts{$split_path} || 0;
            
            if ($expected_reads == $actual_reads) {
                $unchecked->move($split_bam);
                my $meta = $split_bam->metadata;
                $meta->{reads} = $actual_reads;
                $split_bam->add_metadata($meta, replace_data => 1);
                push(@out_files, $split_path);
            }
            else {
                $self->warn("When attempting " . $bam_file->path . " -> " . $split_bam->path . " split, got $actual_reads reads instead of $expected_reads; will delete it");
                $unchecked->remove;
            }
        }
        if (scalar @out_files == scalar keys %split_bams) {
            return 1;
        }
        else {
            $self->throw("Number of correct output files " . scalar @out_files . " not equal to number of expected output files " . scalar keys %split_bams);
        }
    }
}

1;
