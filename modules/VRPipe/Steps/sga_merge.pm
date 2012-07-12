=head1 NAME

VRPipe::Steps::sga_merge - a step

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


# Usage: sga merge [OPTION] ... READS1 READS2
# Merge the sequence files READS1, READS2 into a single file/index
# 
#   -v, --verbose                        display verbose output
#       --help                           display this help and exit
#   -t, --threads=NUM                    use NUM threads to merge the indices (default: 1)
#   -p, --prefix=PREFIX                  write final index to files starting with PREFIX (the default is to concatenate the input filenames)
#   -r, --remove                         remove the original BWT, SAI and reads files after the merge
#   -g, --gap-array=N                    use N bits of storage for each element of the gap array. Acceptable values are 4,8,16 or 32. Lower
#                                        values can substantially reduce the amount of memory required at the cost of less predictable memory usage.
#                                        When this value is set to 32, the memory requirement is essentially deterministic and requires ~5N bytes where
#                                        N is the size of the FM-index of READS2.
#                                        The default value is 4.
#       --no-sequence                    Suppress merging of the sequence files. Use this option when merging the index(es) separate e.g. in parallel
#       --no-forward                     Suppress merging of the forward index. Use this option when merging the index(es) separate e.g. in parallel
#       --no-reverse                     Suppress merging of the reverse index. Use this option when merging the index(es) separate e.g. in parallel

use VRPipe::Base;

class VRPipe::Steps::sga_merge with VRPipe::StepRole {
    method options_definition {
        return { sga_merge_options => VRPipe::StepOption->get(description => 'options to sga index', optional => 1, default_value => ''),
                 sga_exe => VRPipe::StepOption->get(description => 'path to your sga executable', optional => 1, default_value => 'sga') };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => -1, description => 'fastq files to be merged') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            
            my $sga_exe = $options->{sga_exe};
            my $sga_opts = $options->{sga_merge_options};
            if ($sga_opts =~ /merge|-p|--prefix/) {
                $self->throw("sga_merge_options should not include the merge subcommand or --prefix (-p) option");
            }
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga merge '.$sga_opts.' $fastq_file_1 $fastq_file_2'));
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            my @fastqs = sort { $a->s <=> $b->s  } @{$self->inputs->{fastq_files}};
            my $id = 1;
            while (@fastqs) {
                my @fqs;
                # merge the largest file with the smallest file
                push(@fqs, shift @fastqs); 
                push(@fqs, pop @fastqs) if @fastqs;
                my @basenames;
                my $popidx;
                foreach my $fq (@fqs) {
                    my $prefix = $fq->basename;
                    $prefix =~ s/\.(fq|fastq)(\.gz)?//;
                    push @basenames, $prefix;
                    ++$popidx if (-s file($fq->dir, "$prefix.popidx"));
                }
                my $fq_meta = $self->common_metadata(\@fqs);
                my $basename = join '_', @basenames;
                my @outfiles;
                push @outfiles, $self->output_file(output_key => 'merged_fastq_files', basename => "$basename.$id.fq", type => 'fq', metadata => $fq_meta);
                push @outfiles, $self->output_file(output_key => 'merged_bwt_files', basename => "$basename.$id.bwt", type => 'bin', metadata => $fq_meta);
                push @outfiles, $self->output_file(output_key => 'merged_sai_files', basename => "$basename.$id.sai", type => 'txt', metadata => $fq_meta);
                if ($popidx) {
                    push @outfiles, $self->output_file(output_key => 'merged_popidx_files', basename => "$basename.$id.popidx", type => 'txt', metadata => $fq_meta);
                }
                ++$id;
                my $prefix = $outfiles[0]->path;
                $prefix =~ s/\.fq$//;
                my $cmd = qq[$sga_exe merge $sga_opts --prefix $prefix ].join(' ', map { $_->path } @fqs);
                $self->dispatch_wrapped_cmd('VRPipe::Steps::sga_merge', 'merge_and_check', [$cmd, $req, {output_files => \@outfiles}]);
            }
        };
    }
    method outputs_definition {
        return { merged_fastq_files => VRPipe::StepIODefinition->get(type => 'fq', description => 'the files produced by sga index', max_files => -1),
                 merged_bwt_files => VRPipe::StepIODefinition->get(type => 'bin', description => 'the files produced by sga index', max_files => -1),
                 merged_sai_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'the files produced by sga index', max_files => -1),
                 merged_popidx_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'the files produced by sga index', min_files => 0, max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Merge the sequence files READS1, READS2 into a single file/index";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method merge_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($prefix, $in_paths) = $cmd_line =~ /--prefix (\S+) (.+)$/;
        my @in_paths = split ' ', $in_paths;
        @in_paths || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $prefix || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my (@in_fastqs, @in_bwts, @in_sais, @in_popidxs);
        my $popidx;
        foreach my $in_path (@in_paths) {
            push @in_fastqs, VRPipe::File->get(path => $in_path);
            my ($in_bwt, $in_sai, $in_popidx) = ($in_path, $in_path, $in_path);
            $in_bwt =~ s/\.(fq|fastq)/.bwt/;
            push @in_bwts, VRPipe::File->get(path => $in_bwt);
            $in_sai =~ s/\.(fq|fastq)/.sai/;
            push @in_sais, VRPipe::File->get(path => $in_sai);
            $in_popidx =~ s/\.(fq|fastq)/.popidx/;
            if (-s $in_popidx) {
                ++$popidx;
                push @in_popidxs, VRPipe::File->get(path => $in_popidx);
            }
        }
        my $out_fq = VRPipe::File->get(path => $prefix.'.fq');
        my $out_fa = VRPipe::File->get(path => $prefix.'.fa');
        my $out_bwt = VRPipe::File->get(path => $prefix.'.bwt');
        my $out_sai = VRPipe::File->get(path => $prefix.'.sai');
        my $out_popidx;
        if ($popidx) {
            $out_popidx = VRPipe::File->get(path => $prefix.'.popidx');
        }
        
        if (scalar @in_paths == 1) {
            $in_fastqs[0]->copy($out_fq);
            $in_bwts[0]->copy($out_bwt);
            $in_sais[0]->copy($out_sai);
            if ($popidx) {
                $in_popidxs[0]->copy($out_popidx);
            }
            return;
        }
        else {
            $out_fq->disconnect;
            system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        }
        $out_fa->update_stats_from_disc(retries => 3);
        $out_fa->disconnect;
        
        $out_fa->move($out_fq);
        $out_fq->update_stats_from_disc(retries => 3);
        
        my $reads = 0;
        foreach my $fq (@in_fastqs) {
            my $meta = $fq->metadata;
            $reads += $meta->{reads} || $fq->num_records;
        }
        my $actual_reads = $out_fq->num_records;
        
        if ($actual_reads == $reads) {
            $out_fq->add_metadata({reads => $actual_reads});
            return 1;
        }
        else {
            $out_fq->unlink;
            $out_bwt->unlink;
            $out_sai->unlink;
            if ($popidx) {
                $out_popidx->unlink;
            }
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output fastq file, yet there were $reads reads in the original fastq files");
        }
    }
}

1;
