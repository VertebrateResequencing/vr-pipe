
=head1 NAME

VRPipe::Steps::bam_split_by_region - a step

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

class VRPipe::Steps::bam_split_by_region with VRPipe::StepRole {
    use VRPipe::Parser;
    use VRPipe::Parser::bam;
    
    method options_definition {
<<<<<<< Updated upstream
        return {
            reference_index                   => VRPipe::StepOption->create(description => 'Reference index (fai) file.'),
            split_bam_by_region_chunk_size    => VRPipe::StepOption->create(description => 'Chunk size to split the bam into.', default_value => 50000000),
            split_bam_by_region_chunk_overlap => VRPipe::StepOption->create(description => 'Number of bases to allow chunks to overlap.', default_value => 1000),
            split_bam_by_region_chrom_list    => VRPipe::StepOption->create(description => 'Number of bases to allow chunks to overlap.', optional => 1),
            split_bam_by_region_include_mate  => VRPipe::StepOption->create(description => 'boolean; if true mates mapping to the split sequence will be included', optional => 1, default_value => 0)
        };
=======
        return { reference_index                   => VRPipe::StepOption->create(description => 'Reference index (fai) file.'),
                 split_bam_by_region_chunk_size    => VRPipe::StepOption->create(description => 'Chunk size to split the bam into.',                                                    default_value => 50000000),
                 split_bam_by_region_chunk_overlap => VRPipe::StepOption->create(description => 'Number of bases to allow chunks to overlap.',                                          default_value => 1000),
                 split_bam_by_region_chrom_list    => VRPipe::StepOption->create(description => 'Number of bases to allow chunks to overlap.',                           optional => 1),
                 split_bam_by_region_include_mate  => VRPipe::StepOption->create(description => 'boolean; if true mates mapping to the split sequence will be included', optional => 1, default_value => 0) };
>>>>>>> Stashed changes
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to be split') };
    }
    
    method body_sub {
        return sub {
            use VRPipe::Utils::GenomeChunking;
            
            my $self            = shift;
            my $options         = $self->options;
            my $reference_index = $options->{reference_index};
<<<<<<< Updated upstream
            my $chunk_size      = $options->{split_bam_by_region_chunk_size};
            my $chunk_overlap   = $options->{split_bam_by_region_chunk_overlap};
=======
            my $chunk_size      = $options->{split_bam_by_region_chunk_size };
            my $chunk_overlap   = $options->{split_bam_by_region_chunk_overlap };
>>>>>>> Stashed changes
            my $chrom_list      = $options->{chrom_list};
            
            my @chroms = ();
            @chroms = split(' ', $chrom_list) if ($chrom_list);
            
            my $chunk_util = VRPipe::Utils::GenomeChunking->new();
            my $chunks = $chunk_util->chunks(reference_index => $reference_index, chunk_size => $chunk_size, chunk_overlap => $chunk_overlap, chroms => \@chroms);
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $split_dir = $self->output_root->absolute->stringify;
                my $bam_path  = $bam->path->stringify;
                
<<<<<<< Updated upstream
                my $meta = $bam->metadata;
=======
                my $meta      = $bam->metadata;
>>>>>>> Stashed changes
                delete $meta->{reads} if (exists $meta->{reads});
                delete $meta->{bases} if (exists $meta->{bases});
                
                my $chrom_select = $meta->{chrom} || $meta->{split_sequence} || '';
<<<<<<< Updated upstream
                $chrom_select =~ s/^chrom//;
=======
                $chrom_select    =~ s/^chrom//;
>>>>>>> Stashed changes
                
                my (@regions, @output_files);
                my $args = "q[$bam_path], split_dir => q[$split_dir], include_mate => $$options{split_bam_by_region_include_mate}";
                foreach my $chunk (@$chunks) {
                    if ($chrom_select) {
                        my $chrom = $chunk->{chrom};
                        $chrom =~ s/^chrom//;
                        next if ($chrom ne $chrom_select);
                    }
                    my $region = "$$chunk{chrom}_$$chunk{from}-$$chunk{to}";
                    push(@regions, $region);
<<<<<<< Updated upstream
                    push(
                        @output_files,
                        $self->output_file(
                            output_key => 'split_region_bam_files',
                            basename   => "$region." . $bam->basename,
                            type       => 'bam',
                            metadata   => { %$meta, chrom => $$chunk{chrom}, from => $$chunk{from}, to => $$chunk{to} }
                        )
                    );
=======
                    push(@output_files, $self->output_file(output_key => 'split_region_bam_files',
                                                           basename   => "$region." . $bam->basename,
                                                           type       => 'bam',
                                                           metadata   => { %$meta, chrom => $$chunk{chrom}, from => $$chunk{from}, to => $$chunk{to} }));
>>>>>>> Stashed changes
                }
                if (@regions) {
                    my $this_cmd = "use VRPipe::Steps::bam_split_by_region; VRPipe::Steps::bam_split_by_region->split_bam_by_region($args, regions => [qw(@regions)]);";
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [@output_files] });
                }
            }
        };
    }
    
    method outputs_definition {
<<<<<<< Updated upstream
        return {
            split_region_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'split bams for each sequence the reads were aligned to',
                metadata    => {
                    reads => 'number of reads in the split region bam',
                    chrom => 'The chromosome for this split region',
                    from  => 'Split from',
                    to    => 'Split to'
                }
            )
        };
=======
        return { split_region_bam_files => VRPipe::StepIODefinition->create(type        => 'bam',
                                                                            max_files   => -1,
                                                                            description => 'split bams for each sequence the reads were aligned to',
                                                                            metadata    => {
                                                                                             reads => 'number of reads in the split region bam',
                                                                                             chrom => 'The chromosome for this split region',
                                                                                             from  => 'Split from',
                                                                                             to    => 'Split to' }) };
>>>>>>> Stashed changes
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Splits a BAM file into multiple BAM files, one for each sequence the reads were aligned to";
    }
    
    method max_simultaneous {
<<<<<<< Updated upstream
        return 0;          # meaning unlimited
=======
        return 0;            # meaning unlimited
>>>>>>> Stashed changes
    }
    
    method split_bam_by_region (ClassName|Object $self: Str $bam_file!, Str :$split_dir!, ArrayRef :$regions!, Bool :$include_mate = 0) {
        unless (ref($bam_file) && ref($bam_file) eq 'VRPipe::File') {
            $bam_file = VRPipe::File->get(path => file($bam_file));
        }
        my $pb = VRPipe::Parser->create('bam', { file => $bam_file });
        
        my %region_bams;
        foreach my $region (@$regions) {
            my $region_bam = VRPipe::File->create(path => file($split_dir, $region . '.' . $bam_file->basename));
            $region_bams{$region} = $region_bam->path->stringify;
        }
        
        my @get_fields = ('RNAME', 'POS');
        if ($include_mate) {
            push(@get_fields, 'MRNM', 'MPOS');
        }
        
        $bam_file->disconnect;
        
        # stream through the bam, outputting all our desired files
        my $parsed_record = $pb->parsed_record();
        $pb->get_fields(@get_fields);
        my %counts;
        while ($pb->next_record) {
            my $seq_name = $parsed_record->{RNAME};
            my $seq_pos  = $parsed_record->{POS};
            my %out_bams;
            my $seq_bams = $self->position_to_bam($seq_name, $seq_pos, \%region_bams);
            foreach my $seq_bam (@$seq_bams) {
                $out_bams{$seq_bam} = 1;
            }
            if ($include_mate) {
                my $mate_name = $parsed_record->{MRNM};
                my $mate_pos  = $parsed_record->{MPOS};
                my $mate_bams = $self->position_to_bam($mate_name, $mate_pos, \%region_bams);
                foreach my $mate_bam (@$mate_bams) {
                    $out_bams{$mate_bam} = 1;
                }
            }
            
            if (keys %out_bams) {
                foreach my $out_bam (keys %out_bams) {
                    $pb->write_result($out_bam . '.unchecked.bam');
                    $counts{$out_bam}++;
                }
            }
        }
        $pb->close;
        
        # check all the bams we created
        my @out_files;
        foreach my $region_path (values %region_bams) {
            my $unchecked_path = $region_path . '.unchecked.bam';
            my $region_bam     = VRPipe::File->get(path => file($region_path));
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
            my $expected_reads = $counts{$region_path} || 0;
            
            if ($expected_reads == $actual_reads) {
                $unchecked->move($region_bam);
                $region_bam->add_metadata({ reads => $actual_reads }, replace_data => 1);
                push(@out_files, $region_path);
            }
            else {
                $self->warn("When attempting " . $bam_file->path . " -> " . $region_bam->path . " split, got $actual_reads reads instead of $expected_reads; will delete it");
                $unchecked->remove;
            }
        }
        if (scalar @out_files == scalar keys %region_bams) {
            return 1;
        }
        else {
            $self->throw("Number of correct output files " . scalar @out_files . " not equal to number of expected output files " . scalar keys %region_bams);
        }
    }
    
    method position_to_bam (ClassName|Object $self: Str $seq!, Int $pos!, HashRef $region_bams!) {
        my @matching_bams;
        while (my ($region, $bam) = each %$region_bams) {
            my ($chrom, $from, $to) = $region =~ m/^(.+)_(\d+)-(\d+)/;
            next unless ($seq eq $chrom);
            next unless ($pos >= $from);
            next unless ($pos <= $to);
            push @matching_bams, $bam;
        }
        return \@matching_bams;
    }
}

1;
