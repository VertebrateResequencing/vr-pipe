
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
        return {};
    }
    
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->create(type        => 'fq',
                                                                 max_files   => -1,
                                                                 description => 'fastq files to be indexed',
                                                                 metadata    => { sample => 'sample name' }) };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my @fq_paths = map { $_->path } @{ $self->inputs->{fastq_files} };
            my $fq_meta = $self->common_metadata($self->inputs->{fastq_files});
            
            my $merged_fastq_file = $self->output_file(output_key => 'merged_fastq_file', basename => 'merged.fq',     type => 'fq',  metadata => $fq_meta);
            my $index_file        = $self->output_file(output_key => 'index_file',        basename => 'merged.popidx', type => 'txt', metadata => $fq_meta);
            
            my $merged_fastq_path = $merged_fastq_file->path;
            my $index_path        = $index_file->path;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $this_cmd = "use VRPipe::Steps::fastq_merge_and_index; VRPipe::Steps::fastq_merge_and_index->fastq_merge_and_index([qw(@fq_paths)], q[$merged_fastq_path], q[$index_path]);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$merged_fastq_file, $index_file] });
        };
    }
    
    method outputs_definition {
        return { merged_fastq_file => VRPipe::StepIODefinition->create(type => 'fq',  description => 'the files produced by sga index', max_files => 1),
                 index_file        => VRPipe::StepIODefinition->create(type => 'txt', description => 'the files produced by sga index', max_files => 1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merge sample fastq files together and create a popidx file for quick lookup";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method fastq_merge_and_index (ClassName|Object $self: ArrayRef[Str|File] $fastqs!, Str|File $merged_fastq!, Str|File $pop_index!) {
        unless (ref($merged_fastq) && ref($merged_fastq) eq 'VRPipe::File') {
            $merged_fastq = VRPipe::File->get(path => file($merged_fastq));
        }
        unless (ref($pop_index) && ref($pop_index) eq 'VRPipe::File') {
            $pop_index = VRPipe::File->get(path => file($pop_index));
        }
        
        my $seq_fh = $merged_fastq->openw;
        my $idx_fh = $pop_index->openw;
        
        my %samples;
        my $expected_lines = 0;
        foreach my $fq (@$fastqs) {
            my $fq_file = VRPipe::File->get(path => file($fq));
            my $sample = $fq_file->metadata->{sample};
            push @{ $samples{$sample} }, $fq;
            $expected_lines += $fq_file->lines;
        }
        
        $merged_fastq->disconnect;
        
        # Track position data
        my $current_index       = 0;
        my $current_start_index = 0;
        my $current_label       = '';
        
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
                # parse the fastq file
                my $pars = VRPipe::Parser->create('fastq', { file => $fq });
                my $parsed_record = $pars->parsed_record();
                while ($pars->next_record()) {
                    my $id          = $parsed_record->[0];
                    my $seq_string  = $parsed_record->[1];
                    my $qual_string = $parsed_record->[2];
                    
                    my ($header) = split(' ', $id);
                    my $record = "\@$header\n$seq_string\n+\n$qual_string\n";
                    
                    $current_index++;
                    print $seq_fh $record;
                }
            }
        }
        $seq_fh->close;
        
        # Write the last element of the index
        print $idx_fh join("\t", ($current_start_index, $current_index - 1, $current_label)) . "\n";
        $idx_fh->close;
        
        $pop_index->update_stats_from_disc(retries => 3);
        $merged_fastq->update_stats_from_disc(retries => 3);
        
        my $actual_lines = $merged_fastq->lines;
        my $index_lines  = $pop_index->lines;
        my $reads        = $merged_fastq->num_records;
        if ($actual_lines != $expected_lines) {
            $merged_fastq->unlink;
            $pop_index->unlink;
            $self->throw("Merged fastq had $actual_lines actual lines, whereas we expected $expected_lines lines.");
        }
        elsif ($index_lines != scalar @$fastqs) {
            $merged_fastq->unlink;
            $pop_index->unlink;
            $self->throw("Index had $index_lines, whereas we expected " . scalar(@$fastqs) . " lines.");
        }
        else {
            $merged_fastq->add_metadata({ reads => $reads });
            return 1;
        }
    }
}

1;
