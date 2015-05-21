
=head1 NAME

VRPipe::Steps::split_fastq - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk> and Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::split_fastq with VRPipe::StepRole {
    use POSIX;
    
    method options_definition {
        return {
            fastq_chunk_size => VRPipe::StepOption->create(description => 'when splitting fastq files into smaller chunks, this sets the size in bp; a good figure might be 1000000000 for a fast mapper'),
            split_exe        => VRPipe::StepOption->create(
                description   => 'path to your GNU split executable',
                optional      => 1,
                default_value => 'split'
            )
        };
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
            
            my $split_exe  = $self->options->{split_exe};
            my $chunk_size = $self->options->{fastq_chunk_size};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'split',
                    version => VRPipe::StepCmdSummary->determine_version($split_exe . ' --version', '^split \(GNU coreutils\) (.+)$'),
                    summary => "split -l \$lines --filter='gzip -c > \$FILE.gz' \$fastq_file"
                )
            );
            
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
            
            while (my ($ended, $fqs) = each %ended) {
                next unless $fqs;
                my $split_subdir = $ended . '_' . $chunk_size;
                my $split_dir = dir($out_root, $split_subdir);
                
                my @outs = VRPipe::Steps::split_fastq->split_fastq_outs(split_dir => $split_dir, fastqs => $fqs, chunk_size => $chunk_size);
                my @ofiles;
                foreach my $out (@outs) {
                    push(@ofiles, $self->output_file(output_key => 'split_fastq_files', sub_dir => $split_subdir, basename => $out->basename, type => 'fq'));
                }
                
                $self->dispatch_vrpipecode(qq[use VRPipe::Steps::split_fastq; VRPipe::Steps::split_fastq->split_fastq(split_dir => q[$split_dir], split_exe => q[$split_exe], fastqs => [qw(@$fqs)], chunk_size => $chunk_size);], $req, { output_files => \@ofiles });
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

=head2 split_fastq
 
 Title   : split_fastq
 Usage   : $obj->split_fastq(fastqs => \@fastqs,
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
    
    method split_fastq (ClassName|Object $self: Str|Dir :$split_dir!, Str :$split_exe!, ArrayRef[Str|File] :$fastqs!, PositiveInt :$chunk_size!) {
        -d $split_dir || $self->throw("split_fastq split_dir '$split_dir' does not exist");
        
        my $splits = $self->_splits($fastqs, $chunk_size);
        
        my @fqs;
        my $read_lengths = 0;
        my %parent_metadata;
        foreach my $fq_path (@$fastqs) {
            my $basename = file($fq_path)->basename;
            my $prefix   = $basename;
            $prefix =~ s/\.f[^.]+(?:\.gz)?$//;
            
            my $suffix = 'fastq.gz';
            if ($splits == 1 && $basename !~ /\.gz$/) {
                $suffix = 'fastq';
            }
            
            my $fq_file = VRPipe::File->get(path => $fq_path);
            $fq_file->update_stats_from_disc(); # necessary in case files are deleted and then we redo this step
            my $fq_meta = $fq_file->metadata;
            $read_lengths += $fq_meta->{avg_read_length};
            
            my $source_path = file($fq_path)->resolve->stringify;
            $parent_metadata{$source_path} = $fq_meta;
            push(@fqs, [$prefix, $fq_file, $suffix]);
        
        }
        
        # just symlink if we only have 1 chunk
        if ($splits == 1) {
            foreach my $i (0 .. $#fqs) {
                my ($prefix, $fq_file, $suffix) = @{ $fqs[$i] };
                my $split_file = VRPipe::File->get(path => file($split_dir, "$prefix.000.$suffix")->stringify, type => 'fq');
                symlink($fq_file->path, $split_file->path); # don't use $fq_file->symlink because $split_file->metadata must be different
                $split_file->update_stats_from_disc;
                my $meta = $fq_file->metadata;
                delete $meta->{expected_md5};
                $meta->{chunk} = 1;
                $split_file->add_metadata($meta);
            }
            return 1;
        }
        
        my $seqs_per_split = floor($chunk_size / $read_lengths);
        
        $fqs[0]->[1]->disconnect;
        foreach my $i (0 .. $#fqs) {
            my ($prefix, $fq_file, $suffix) = @{ $fqs[$i] };
            my $fq_path = $fq_file->path;
            
            my $cmd_line;
            my $lines_per_split = 4 * $seqs_per_split;
            if ($fq_path =~ /\.gz$/) {
                $cmd_line = "gunzip -c $fq_path | $split_exe -da 3 -l $lines_per_split --filter='gzip -c > \$FILE.gz' - --additional-suffix=.fastq $prefix.";
            }
            else {
                $cmd_line = "cat $fq_path | $split_exe -da 3 -l $lines_per_split --filter='gzip -c > \$FILE.gz' - --additional-suffix=.fastq $prefix.";
            }
            system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        }
        
        my %seq_base_counts;
        # check the chunks seem fine and construct metadata
        foreach my $i (0 .. $#fqs) {
            my ($prefix, $fq_file, $suffix) = @{ $fqs[$i] };
            my $in_lines    = $fq_file->lines;
            my $fq_meta     = $fq_file->metadata;
            my $read_length = $fq_meta->{avg_read_length};
            
            my $out_lines = 0;
            my @split_files;
            foreach my $split_num ("000" .. sprintf("%03d", $splits - 1)) {
                my $split_file = VRPipe::File->get(path => file($split_dir, "$prefix." . sprintf("%03d", $split_num) . ".$suffix")->stringify, type => 'fq');
                $split_file->update_stats_from_disc();
                my $seqs        = $split_file->lines;
                my $base_counts = $seqs * $read_length;
                $out_lines += $seqs;
                my $mate;
                
                if (@fqs == 2) {
                    my $mate_ref;
                    if ($i == 0) {
                        $mate_ref = $fqs[1];
                    }
                    else {
                        $mate_ref = $fqs[0];
                    }
                    $mate = file($split_dir, "$mate_ref->[0]." . sprintf("%03d", $split_num) . ".$suffix")->stringify;
                }
                
                $seq_base_counts{ file($split_dir, "$prefix." . sprintf("%03d", $split_num) . ".$suffix")->stringify } = [$seqs, $base_counts, $mate, $split_num];
                push(@split_files, $split_file->path);
            }
            
            unless ($out_lines == $in_lines) {
                $self->throw($fq_file->path . " had $in_lines lines, but the split files @split_files ended up with only $out_lines!");
            }
        }
        
        # add metadata to each output file
        while (my ($path, $stats) = each %seq_base_counts) {
            my ($seqs, $num_bases, $mate, $chunk) = @{$stats};
            my $avg_read_length = sprintf("%0.2f", $num_bases / $seqs);
            
            my $split_file  = VRPipe::File->get(path => $path);
            my $split_meta  = $split_file->metadata;
            my $parent_meta = $parent_metadata{ $split_meta->{source_fastq} };
            $self->throw("no parent metadata for parent " . $split_meta->{source_fastq}) unless ($parent_meta && defined $parent_meta->{paired});
            
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
        
        return 1;
    }
    
    method split_fastq_outs (ClassName|Object $self: Str|Dir :$split_dir!, ArrayRef[Str|File] :$fastqs!, PositiveInt :$chunk_size!) {
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
                        path => file($split_dir, "$prefix.000.$suffix"),
                        type => 'fq',
                        metadata => { source_fastq => $source_fastq }
                    )
                );
            }
            else {
                for my $split_num ("000" .. sprintf("%03d", $splits - 1)) {
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
