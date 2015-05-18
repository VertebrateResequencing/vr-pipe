
=head1 NAME

VRPipe::Steps::samtools_merge - a step

=head1 DESCRIPTION

Runs the samtools merge.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::samtools_merge with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe           => VRPipe::StepOption->create(description => 'path to samtools executable', optional => 1, default_value => 'samtools'),
            samtools_merge_options => VRPipe::StepOption->create(description => 'options for samtools merge',  optional => 1, default_value => ''),
        };
    }
    
    method inputs_definition {
        return {
            aln_files => VRPipe::StepIODefinition->create(
                type        => 'aln',                                          # cram or bam
                max_files   => -1,
                description => '1 or more coordinate sorted BAM or CRAM files',
                metadata    => {
                    bases    => 'total number of base pairs',
                    reads    => 'total number of reads (sequences)',
                    paired   => '0=single ended reads only; 1=paired end reads present',
                    optional => ['bases', 'reads', 'paired']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                = shift;
            my $options             = $self->options;
            my $samtools            = $options->{samtools_exe};
            my $samtools_merge_opts = $options->{samtools_merge_options};
            my $check_read_sum      = $samtools_merge_opts =~ m/-R/ ? 0 : 1;
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools merge $samtools_merge_opts \$merged_bam \$input_files(s); samtools index \$merged_bam"
                )
            );
            
            my $req = $self->new_requirements(memory => 3000, time => 1);
            
            my $inputs          = $self->inputs->{aln_files};
            my $merged_metadata = $self->common_metadata($inputs);
            $merged_metadata = { %$merged_metadata, $self->element_meta };
            my @input_paths = map { $_->path } @$inputs;
            my $file_list = VRPipe::FileList->create(files => $inputs)->id;
            my $basename = $self->step_state->dataelement->id;
            if (defined $$merged_metadata{chrom} && defined $$merged_metadata{from} && defined $$merged_metadata{to}) {
                my ($chrom, $from, $to) = ($$merged_metadata{chrom}, $$merged_metadata{from}, $$merged_metadata{to});
                $samtools_merge_opts .= " -R $chrom:$from-$to";
                $basename       = "${chrom}_${from}-${to}.$basename";
                $check_read_sum = 0;
            }
            my $merged_bam_file = $self->output_file(
                output_key => 'merged_bam_files',
                basename   => $basename . '.bam',
                type       => 'bam',
                metadata   => $merged_metadata
            );
            my $merged_bai_file = $self->output_file(
                output_key => 'merged_bai_files',
                basename   => $basename . '.bam.bai',
                type       => 'bai',
                metadata   => $merged_metadata
            );
            my $merged_path = $merged_bam_file->path;
            my $merge_cmd   = qq[$samtools merge $samtools_merge_opts $merged_path @input_paths];
            my $index_cmd   = qq[$samtools index $merged_path];
            my $this_cmd    = "use VRPipe::Steps::samtools_merge; VRPipe::Steps::samtools_merge->merge_and_check(qq[$merge_cmd; $index_cmd], input_file_list => $file_list, out_file => q[$merged_path], check_read_sum => $check_read_sum);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$merged_bam_file, $merged_bai_file] });
        };
    }
    
    method outputs_definition {
        return {
            merged_bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => 1, description => 'a merged BAM file'),
            merged_bai_files => VRPipe::StepIODefinition->create(type => 'bai', max_files => 1, description => 'BAI index file for merged BAM file'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges CRAM and BAM files using samtools merge to create a merged BAM file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method merge_and_check (ClassName|Object $self: Str $cmd_line, Int :$input_file_list!, Str|File :$out_file!, Bool :$check_read_sum = 1 ) {
        unless (ref($out_file) && ref($out_file) eq 'VRPipe::File') {
            $out_file = VRPipe::File->create(path => file($out_file));
        }
        
        my $out_path  = $out_file->path;
        my $out_index = VRPipe::File->get(path => $out_file->type eq 'bam' ? qq[$out_path.bai] : qq[$out_path.crai]);
        my @in_files  = VRPipe::FileList->get(id => $input_file_list)->files;
        
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        $out_index->update_stats_from_disc(retries => 3);
        
        my $actual_reads = $out_file->num_records;
        if ($check_read_sum) {
            my $expected_reads  = 0;
            my $expected_bases  = 0;
            my $expected_paired = 0;
            my $paired          = 0;
            foreach my $in_file (@in_files) {
                my $meta = $in_file->metadata;
                $in_file->disconnect;
                $expected_reads += $meta->{reads} || $in_file->num_records;
                $expected_bases += $meta->{bases} if $meta->{bases};
                $paired = 1 if (exists $meta->{paired});
                $expected_paired ||= $meta->{paired};
            }
            
            if ($actual_reads == $expected_reads) {
                my %new_meta;
                $new_meta{reads}  = $actual_reads;
                $new_meta{bases}  = $expected_bases if $expected_bases;
                $new_meta{paired} = $expected_paired if $paired;
                $out_file->add_metadata(\%new_meta);
                return 1;
            }
            else {
                $out_file->unlink;
                $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output CRAM file, yet there were $expected_reads reads in the original CRAM file");
            }
        }
        else {
            return 1;
        }
    }
}

1;
