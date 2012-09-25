
=head1 NAME

VRPipe::Steps::fastq_merge_and_index - a step

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

class VRPipe::Steps::fastq_merge_and_index with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { fastq_merge_and_index_compress_fastq => VRPipe::StepOption->create(description => 'compress the fastq output of fastq_merge_and_index (boolean)', optional => 1, default_value => 1) };
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'fastq files to be merged',
                metadata    => {
                    sample => 'sample name',
                    reads  => 'number of reads in the fastq',
                    bases  => 'number of bases in the fastq'
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $compress = $options->{fastq_merge_and_index_compress_fastq};
            
            my $fq_meta = $self->common_metadata($self->inputs->{fastq_files});
            
            my $merge_list = $self->output_file(basename => "merge_list.txt", type => 'txt', temporary => 1);
            my $merge_list_path = $merge_list->path;
            $merge_list->create_fofn($self->inputs->{fastq_files});
            
            my ($reads, $bases) = (0, 0);
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                my $meta = $fq->metadata;
                $reads += $meta->{reads};
                $bases += $meta->{bases};
            }
            $fq_meta->{reads} = $reads;
            $fq_meta->{bases} = $bases;
            
            my $fq_basename = 'merged.fq';
            if ($compress) {
                $fq_basename .= '.gz';
            }
            my $merged_fastq_file = $self->output_file(output_key => 'merged_fastq_file', basename => $fq_basename,    type => 'fq',  metadata => $fq_meta);
            my $index_file        = $self->output_file(output_key => 'index_file',        basename => 'merged.popidx', type => 'txt', metadata => $fq_meta);
            
            my $merged_fastq_path = $merged_fastq_file->path;
            my $index_path        = $index_file->path;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $this_cmd = "use VRPipe::Steps::fastq_merge_and_index; VRPipe::Steps::fastq_merge_and_index->fastq_merge_and_index(fastq_fofn => q[$merge_list_path], merged_fastq => q[$merged_fastq_path], pop_index => q[$index_path]);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$merged_fastq_file, $index_file] });
        };
    }
    
    method outputs_definition {
        return {
            merged_fastq_file => VRPipe::StepIODefinition->create(
                type        => 'fq',
                description => 'the merged fastq file',
                max_files   => 1,
                metadata    => {
                    reads => 'number of reads in the merged fastq',
                    bases => 'number of bases in the merged fastq'
                }
            ),
            index_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'popidx file to map sequences back to samples', max_files => 1)
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merge sample fastq files together and create a popidx file for quick lookup";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method fastq_merge_and_index (ClassName|Object $self: Str|File :$fastq_fofn!, Str|File :$merged_fastq!, Str|File :$pop_index!) {
        unless (ref($fastq_fofn) && ref($fastq_fofn) eq 'VRPipe::File') {
            $fastq_fofn = VRPipe::File->get(path => file($fastq_fofn));
        }
        unless (ref($merged_fastq) && ref($merged_fastq) eq 'VRPipe::File') {
            $merged_fastq = VRPipe::File->get(path => file($merged_fastq));
        }
        unless (ref($pop_index) && ref($pop_index) eq 'VRPipe::File') {
            $pop_index = VRPipe::File->get(path => file($pop_index));
        }
        
        my $seq_fh = $merged_fastq->openw;
        my $idx_fh = $pop_index->openw;
        
        my (%samples, %reads);
        my $expected_lines = 0;
        my $fq_fh          = $fastq_fofn->openr;
        while (my $fq_path = <$fq_fh>) {
            chomp $fq_path;
            my $fq_file = VRPipe::File->get(path => $fq_path);
            my $sample = $fq_file->metadata->{sample};
            push @{ $samples{$sample} }, $fq_path;
            $reads{$fq_path} = $fq_file->metadata->{reads};
            $expected_lines += $fq_file->lines(raw => 1);
        }
        $fq_fh->close;
        
        # Track position data
        my $current_index       = 0;
        my $current_start_index = 0;
        my $current_label       = '';
        my $written_lines       = 0;
        
        $merged_fastq->disconnect;
        
        # Iterate over every file
        while (my ($sample, $fqs) = each %samples) {
            # Write out the previous file's data to the index, if any
            if ($current_label ne '') {
                print $idx_fh join("\t", ($current_start_index, $current_index - 1, $current_label)) . "\n";
            }
            
            # Reset position data
            $current_label       = $sample;
            $current_start_index = $current_index;
            
            foreach my $fq (@$fqs) {
                # write the contents of fq to the merged fq file
                # update the position of the current index
                my $fh;
                if ($fq =~ /\.gz$/) {
                    open($fh, "gunzip -c $fq |") or $self->throw("Could not open 'gunzip -c $fq |': $!");
                }
                else {
                    open($fh, "< $fq") or $self->throw("Could not open '< $fq': $!");
                }
                $fh || $self->throw("Could get filehandle for file $fq");
                while (my $line = <$fh>) {
                    print $seq_fh $line;
                    $written_lines++;
                }
                close $fh;
                $current_index += $reads{$fq};
            }
        }
        $seq_fh->close;
        
        # Write the last element of the index
        print $idx_fh join("\t", ($current_start_index, $current_index - 1, $current_label)) . "\n";
        $idx_fh->close;
        
        $pop_index->update_stats_from_disc(retries => 3);
        $merged_fastq->update_stats_from_disc(retries => 3);
        
        my $index_lines = $pop_index->lines;
        if ($written_lines != $expected_lines) {
            $merged_fastq->unlink;
            $pop_index->unlink;
            $self->throw("Merged fastq had $written_lines written lines, whereas we expected $expected_lines lines.");
        }
        elsif ($index_lines != $fastq_fofn->lines) {
            $merged_fastq->unlink;
            $pop_index->unlink;
            $self->throw("Index had $index_lines, whereas we expected " . $fastq_fofn->lines . " lines.");
        }
        else {
            return 1;
        }
    }
}

1;
