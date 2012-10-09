
=head1 NAME

VRPipe::Steps::fastq_split - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::fastq_split with VRPipe::StepRole {
    use POSIX;
    
    method options_definition {
        return { fastq_chunk_size => VRPipe::StepOption->create(description => 'when splitting fastq files into smaller chunks, this sets the size in bp; a good figure might be 1000000000 for a fast mapper') };
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => 3,
                description => '1-3 fastq files of the same lane',
                metadata    => {
                    lane            => 'lane name (a unique identifer for this sequencing run)',
                    bases           => 'total number of base pairs',
                    reads           => 'total number of reads (sequences)',
                    avg_read_length => 'the average length of reads',
                    paired          => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                    mate            => 'if paired, the path to the fastq that is our mate',
                    optional        => ['mate']
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my ($se, @pe, %ended);
            foreach my $fastq (@{ $self->inputs->{fastq_files} }) {
                my $metadata = $fastq->metadata;
                my $paired   = $metadata->{paired};
                if ($paired == 0) {
                    $se = $fastq->path;
                }
                elsif ($paired == 1) {
                    unshift(@pe, $fastq->path);
                }
                elsif ($paired == 2) {
                    push(@pe, $fastq->path);
                }
            }
            $ended{pe} = @pe == 2 ? \@pe  : undef;
            $ended{se} = $se      ? [$se] : undef;
            my $out_root = $self->output_root;
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my $chunk_size = $self->options->{fastq_chunk_size};
            
            while (my ($ended, $fqs) = each %ended) {
                next unless $fqs;
                my $split_subdir = $ended . '_' . $chunk_size;
                my $split_dir = Path::Class::Dir->new($out_root, $split_subdir);
                
                my @outs = VRPipe::Steps::fastq_split->fastq_split_outs(split_dir => $split_dir, fastqs => $fqs, chunk_size => $chunk_size);
                foreach my $out (@outs) {
                    $self->output_file(output_key => 'split_fastq_files', sub_dir => $split_subdir, basename => $out->basename, type => 'fq');
                }
                
                $self->dispatch_vrpipecode(qq[use VRPipe::Steps::fastq_split; VRPipe::Steps::fastq_split->fastq_split(split_dir => q[$split_dir], fastqs => [qw(@$fqs)], chunk_size => $chunk_size);], $req);
            }
        };
    }
    
    method outputs_definition {
        return {
            split_fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'split fastq files',
                metadata    => {
                    source_fastq    => 'the fastq file this was split from',
                    bases           => 'total number of base pairs',
                    reads           => 'total number of reads (sequences) in this chunk',
                    avg_read_length => 'the average length of reads',
                    paired          => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                    mate            => 'if paired, the path to the fastq that is our mate and of the corresponding chunk',
                    chunk           => 'an int to say which chunk of the source_fastq this is',
                    optional        => ['mate']
                }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Takes a single-ended fastq file and/or 2 fastq files of a paired-end run and splits them into multiple smaller fastq files";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }

=head2 fastq_split
 
 Title   : fastq_split
 Usage   : $obj->fastq_split(fastqs => \@fastqs,
                             split_dir => '/path/to/desired/split_dir',
                             chunk_size => 1000000);
 Function: Split the fastq(s) into multiple smaller files. If more than one
           fastq file is supplied, they are treated as a unit with regard to
           the chunk size. So if you had two fastq files where all sequences
           were 100bp long, and you supplied only one of them with a chunk_size
           of 1000, you'd end up with chunks containing 10 sequences each. But
           with both fastq files supplied, you'd end up with chunks containing
           5 sequences each. The idea being, you'll then use the Nth chunk of
           both fastqs at the same time, bringing the total bases to the
           chunk size.
 Returns : int (the number of splits created)
 Args    : fastqs => array ref of fastq files,
           split_dir => '/path/to/desired/split_dir' (location to store the
                                                      resulting files)
           chunk_size => int (max number of bases per chunk)

=cut
    
    method fastq_split (ClassName|Object $self: Str|Dir :$split_dir!, ArrayRef[Str|File] :$fastqs!, PositiveInt :$chunk_size!) {
        -d $split_dir || $self->throw("fastq_split split_dir '$split_dir' does not exist");
        
        my $splits = $self->_splits($fastqs, $chunk_size);
        
        my @ins;
        my @outs;
        my $split_num = 1;
        my %parent_metadata;
        foreach my $fq_path (@$fastqs) {
            my $basename = file($fq_path)->basename;
            my $prefix   = $basename;
            $prefix =~ s/\.f[^.]+(?:\.gz)?$//;
            
            my $fq_file = VRPipe::File->get(path => $fq_path);
            $fq_file->update_stats_from_disc(); # necessary in case files are deleted and then we redo this step
            $parent_metadata{$fq_path} = $fq_file->metadata;
            my $ifh = $fq_file->openr;
            push(@ins, [$ifh, $fq_file]);
            
            my $suffix = 'fastq.gz';
            if ($splits == 1 && $basename !~ /\.gz$/) {
                $suffix = 'fastq';
            }
            
            my $split_file = VRPipe::File->get(path => file($split_dir, "$prefix.$split_num.$suffix"), type => 'fq');
            my $ofh = $split_file->openw;
            push(@outs, [$prefix, $split_file, $ofh]);
        }
        
        # just symlink if we only have 1 chunk
        if ($splits == 1) {
            foreach my $i (0 .. $#ins) {
                my $fq_file = $ins[$i]->[1];
                $fq_file->close();
                my $split_file = $outs[$i]->[1];
                $split_file->close();
                $split_file->unlink;
                symlink($fq_file->path, $split_file->path); # don't use $fq_file->symlink because $split_file->metadata must be different
                
                my $meta = $fq_file->metadata;
                delete $meta->{expected_md5};
                $meta->{chunk} = 1;
                $split_file->add_metadata($meta);
            }
            return 1;
        }
        
        my $num_bases      = 0;
        my $expected_lines = @ins * 4;
        my %base_counts;
        my %seq_base_counts;
        my $seqs = 0;
        while (1) {
            # get the next entry (4 lines) from each input fastq
            my @seqs;
            my $lines       = 0;
            my $these_bases = 0;
            my %these_base_counts;
            foreach my $i (0 .. $#ins) {
                my $in_fh = $ins[$i]->[0];
                
                for (1 .. 4) {
                    my $line = <$in_fh>;
                    defined $line || next;
                    $lines++;
                    
                    push(@{ $seqs[$i] }, $line);
                    
                    if ($_ == 2) {
                        my $seq = $seqs[$i]->[1];
                        chomp($seq);
                        my $seq_length = length($seq);
                        $these_bases += $seq_length;
                        $these_base_counts{$i} += $seq_length;
                    }
                }
            }
            
            # check for truncation/ eof
            if ($lines == 0) {
                last;
            }
            elsif ($lines != $expected_lines) {
                $self->throw("one of the fastq files ended early");
            }
            
            $seqs++;
            
            # start a new chunk if necessary
            $num_bases += $these_bases;
            if ($num_bases > $chunk_size) {
                $split_num++;
                
                foreach my $i (0 .. $#outs) {
                    my $ref = $outs[$i];
                    my ($prefix, $old) = @{$ref};
                    $old->close;
                    
                    my $mate;
                    my $old_split = $split_num - 1;
                    if (@outs == 2) {
                        my $mate_ref;
                        if ($i == 0) {
                            $mate_ref = $outs[1];
                        }
                        else {
                            $mate_ref = $outs[0];
                        }
                        
                        $mate = file($split_dir, "$mate_ref->[0].$old_split.fastq.gz")->stringify;
                    }
                    $seq_base_counts{ file($split_dir, "$prefix.$old_split.fastq.gz")->stringify } = [$seqs - 1, $base_counts{$i}, $mate, $old_split];
                    
                    my $split_file = VRPipe::File->get(path => file($split_dir, "$prefix.$split_num.fastq.gz"), type => 'fq');
                    $ref->[1]        = $split_file;
                    $ref->[2]        = $split_file->openw;
                    $base_counts{$i} = 0;
                }
                
                $num_bases = $these_bases;
                $seqs      = 1;
            }
            
            # print out the entries
            foreach my $i (0 .. $#seqs) {
                my @lines  = @{ $seqs[$i] };
                my $out_fh = $outs[$i]->[2];
                foreach (@lines) {
                    print $out_fh $_;
                }
                
                $base_counts{$i} += $these_base_counts{$i};
            }
        }
        foreach my $ref (@ins, @outs) {
            $ref->[1]->close;
        }
        
        if ($seqs) {
            foreach my $i (0 .. $#outs) {
                my $ref = $outs[$i];
                my ($prefix, $old) = @{$ref};
                $old->close;
                
                my $mate;
                if (@outs == 2) {
                    my $mate_ref;
                    if ($i == 0) {
                        $mate_ref = $outs[1];
                    }
                    else {
                        $mate_ref = $outs[0];
                    }
                    
                    $mate = file($split_dir, "$mate_ref->[0].$split_num.fastq.gz")->stringify;
                }
                $seq_base_counts{ file($split_dir, "$prefix.$split_num.fastq.gz")->stringify } = [$seqs, $base_counts{$i}, $mate, $split_num];
            }
        }
        
        # check the chunks seem fine
        foreach my $i (0 .. $#ins) {
            my $fastq_file = $ins[$i]->[1];
            my $in_lines   = $fastq_file->lines;
            
            my $prefix    = $outs[$i]->[0];
            my $out_lines = 0;
            foreach my $test_split_num (1 .. $split_num) {
                my $split_file = VRPipe::File->get(path => file($split_dir, "$prefix.$test_split_num.fastq.gz"));
                $out_lines += $split_file->lines;
            }
            
            unless ($out_lines == $in_lines) {
                $self->throw("$fastq_file had $in_lines lines, but the split files ended up with only $out_lines!");
            }
        }
        
        unless ($split_num == $splits) {
            $self->throw("Something weird happened: made $split_num splits, but expected $splits!");
        }
        
        # add metadata to each output file
        while (my ($path, $stats) = each %seq_base_counts) {
            my ($seqs, $num_bases, $mate, $chunk) = @{$stats};
            my $avg_read_length = sprintf("%0.2f", $num_bases / $seqs);
            
            my $split_file  = VRPipe::File->get(path => $path);
            my $split_meta  = $split_file->metadata;
            my $parent_meta = $parent_metadata{ $split_meta->{source_fastq} };
            
            while (my ($key, $val) = each %$parent_meta) {
                $split_meta->{$key} = $val;
            }
            
            my $paired = $parent_meta->{paired};
            if ($parent_meta->{mate}) {
                if ($paired == 0) {
                    delete $split_meta->{mate};
                }
                else {
                    $split_meta->{mate} = $mate;
                }
            }
            
            $split_meta->{bases}           = $num_bases;
            $split_meta->{reads}           = $seqs;
            $split_meta->{avg_read_length} = $avg_read_length;
            delete $split_meta->{expected_md5};
            $split_meta->{chunk} = $chunk;
            
            $split_file->add_metadata($split_meta);
        }
        
        return $split_num;
    }
    
    method fastq_split_outs (ClassName|Object $self: Str|Dir :$split_dir!, ArrayRef[Str|File] :$fastqs!, PositiveInt :$chunk_size!) {
        my $splits = $self->_splits($fastqs, $chunk_size);
        
        my @outs;
        foreach my $fq (@$fastqs) {
            my $basename = $fq->basename;
            my $prefix   = $basename;
            $prefix =~ s/\.f[^.]+(?:\.gz)?$//;
            
            # we must not $fq->resolve->stringify because this actually
            # changes $fq, so subsequent things look at $fastqs may encounter
            # files not in the VRPipe db (and so with no metadata)
            my $source_fastq = file($fq)->resolve->stringify;
            
            if ($splits == 1) {
                my $suffix = $fq =~ /\.gz$/ ? 'fastq.gz' : 'fastq';
                push(
                    @outs,
                    VRPipe::File->create(
                        path => file($split_dir, "$prefix.1.$suffix"),
                        type => 'fq',
                        metadata => { source_fastq => $source_fastq }
                    )
                );
            }
            else {
                for my $split_num (1 .. $splits) {
                    push(
                        @outs,
                        VRPipe::File->create(
                            path => file($split_dir, "$prefix.$split_num.fastq.gz"),
                            type => 'fq',
                            metadata => { source_fastq => $source_fastq }
                        )
                    );
                }
            }
        }
        
        return @outs;
    }
    
    method _splits (ClassName|Object $self: ArrayRef[Str|File] $fastqs, PositiveInt $chunk_size) {
        my $read_lengths = 0;
        my $seqs         = 0;
        foreach my $fq (@$fastqs) {
            my $vrfile_meta = VRPipe::File->get(path => $fq)->metadata;
            $read_lengths += $vrfile_meta->{avg_read_length};
            $seqs = $vrfile_meta->{reads};
        }
        my $seqs_per_split = floor($chunk_size / $read_lengths);
        my $splits         = ceil($seqs / $seqs_per_split);
        
        return $splits;
    }
}

1;
