
=head1 NAME

VRPipe::Steps::bam_to_fastq - a step

=head1 DESCRIPTION

Runs the bam2fastq executable to extract sequences from a BAM file into fastq

=head1 AUTHOR

Chris Joyce <ch5@sanger.ac.uk>.

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

class VRPipe::Steps::bam_to_fastq with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { 
            bam2fastq_exe => VRPipe::StepOption->create(description => 'path to bam2fastq executable', optional => 1, default_value => 'bam2fastq'),
            bam2fastq_opts => VRPipe::StepOption->create(description => 'bam2fastq options excluding --o', optional => 1),
            fastqcheck_exe => VRPipe::StepOption->create(description => 'path to fastqcheck executable', optional => 1, default_value => 'fastqcheck'),
        }
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more BAM files',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    forward_reads    => 'number of forward reads',
                    reverse_reads    => 'number of reverse reads',
                    paired           => '0=single ended reads only; 1=paired end reads present',
                    mean_insert_size => 'mean insert size (0 if unpaired)',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study, put in the DS field of the RG header line',
                    optional         => ['library', 'sample', 'center_name', 'platform', 'study', 'mean_insert_size']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                 = shift;
            my $options              = $self->options;
            my $bam2fastq_exe = $options->{bam2fastq_exe};
            my $bam2fastq_opts = $options->{bam2fastq_opts};
            my $fastqcheck_exe = $options->{fastqcheck_exe};
            $bam2fastq_opts  .= ' -f ';   # force overwrite
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta = $bam->metadata;
                next unless $meta->{reads};
                my $paired = $meta->{paired};
                
                my $source_bam = $bam->path->stringify;
                my $fastq_meta = { source_bam => $source_bam };
                foreach my $key (qw(lane insert_size mean_insert_size library sample center_name platform study split_sequence population continent chrom from to)) {
                    if (defined $meta->{$key}) {
                        $fastq_meta->{$key} = $meta->{$key};
                    }
                }
                my $basename = $bam->basename;
                $basename =~ s/\.bam//;
                my $out_spec;
                my @outfiles;

                if ($paired) {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename => "${basename}_1.fastq",
                        type => 'fq',
                        metadata => {
                            %$fastq_meta,
                            reads => $meta->{forward_reads},
                            paired => 1
                        }
                    );
                    my $reverse = $self->output_file(
                        output_key => 'fastq_files',
                        basename => "${basename}_2.fastq",
                        type => 'fq',
                        metadata => {
                            %$fastq_meta,
                            reads => $meta->{reverse_reads},
                            paired => 2
                        }
                    );
                    @outfiles = ($fastq, $reverse);
                    
                    $fastq->add_metadata({ mate => $reverse->path->stringify });
                    $reverse->add_metadata({ mate => $fastq->path->stringify });
                    
                    $out_spec = 'forward => q[' . $fastq->path . '], reverse => q[' . $reverse->path . ']';
                }
                else {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename => "${basename}_M.fastq",
                        type => 'fq',
                        metadata => {
                            %$fastq_meta,
                            reads => $meta->{reads},
                            bases => $meta->{bases},
                            avg_read_length => sprintf("%0.2f", $meta->{reads} / $meta->{bases}),
                            paired => 0
                        }
                    );
                    @outfiles = ($fastq);
                    $out_spec = 'single => q[' . $fastq->path . ']';
                }
                
                my $out_log = $self->output_file(output_key => 'bam2fastq_logs', basename => "$basename.log", type => 'txt');
                push (@outfiles, $out_log);

                my $this_cmd = "use VRPipe::Steps::bam_to_fastq; VRPipe::Steps::bam_to_fastq->bam_to_fastq(bam => q[$source_bam], $out_spec, bam2fastq_exe => q[$bam2fastq_exe], bam2fastq_opts => q[$bam2fastq_opts], fastqcheck_exe => q[$fastqcheck_exe]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
            }
        };
    }
    
    method outputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => '1 or more fastq files',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    avg_read_length  => 'the average length of reads',
                    paired           => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                    mate             => 'if paired, the path to the fastq that is our mate',
                    source_bam       => 'path of the bam file this fastq file was made from',
                    insert_size      => 'expected library insert size (0 if unpaired)',
                    mean_insert_size => 'calculated mean insert size (0 if unpaired)',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study',
                    optional         => ['mate', 'library', 'sample', 'center_name', 'platform', 'study', 'insert_size', 'mean_insert_size']
                }
            ),
            bam2fastq_logs => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'bam2fastq log',
                max_files   => -1
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts bam files to fastq files";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method bam_to_fastq (ClassName|Object $self: Str|File :$bam!, Str|File :$forward?, Str|File :$reverse?, Str|File :$single?, Str|File :$bam2fastq_exe, Str :$bam2fastq_opts, Str|File :$fastqcheck_exe) {

        my $in_file = VRPipe::File->get(path => $bam);
        my $bam_meta = $in_file->metadata;

        # get the output path so we can specify the logfile name and -o param
        my $out_dir;
        if ( $bam_meta->{paired} ) {
            my $fq = VRPipe::File->get(path => $forward);
            $out_dir = $fq->dir; 
        }
        else {
            my $fq = VRPipe::File->get(path => $single);
            $out_dir = $fq->dir unless $out_dir; 
        }
        
        my $basename = $in_file->basename;
        $basename =~ s/\.bam//;
        my $logfile = "$out_dir/$basename.log";
        # '#' is a special char, (replaced with _1 and _2 for PE reads, _M for SE reads).
        my $out_param = "$out_dir/${basename}#.fastq";

        my $cmd = "$bam2fastq_exe $bam2fastq_opts -o $out_param $bam 2>$logfile";
        system($cmd) && $self->throw("failed to run [$cmd]");

        # determine which output files we need to check
        my @out_files;
        if ( $bam_meta->{paired} ) {
            push(@out_files,  VRPipe::File->get(path => $forward));
            push(@out_files,  VRPipe::File->get(path => $reverse));
        }
        else {
            push(@out_files,  VRPipe::File->get(path => $single));
        }

        # Create additional files fastq if necessary, or get rid of the unneeded empty files
        if ( $bam_meta->{paired} ) {
            my $path = $forward;
            $path =~ s/_1.fastq/_M.fastq/;
            add_or_rm_fastq($self,$path,$forward,\@out_files);

        }
        else {
            foreach my $n (1,2) {
                my $path = $single;
                $path =~ s/_M.fastq/_$n.fastq/;
                add_or_rm_fastq($self,$path,$single,\@out_files);
            }
        }

        # check logfile read count agrees with bam meta
        my $total_reads_parsed;
        my $log = VRPipe::File->get(path => "$logfile");
        $log->update_stats_from_disc;
        my $logh = $log->openr;
        while (my $line = <$logh>) {
            $total_reads_parsed = $1 if $line =~ /(\S+) sequences in the BAM file/;
        }
        unless ($total_reads_parsed) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("Failed to get sequence count from $logfile, or count is zero");
        }
        unless ($total_reads_parsed == $bam_meta->{reads}) {
            foreach my $out_file (@out_files) {
                $out_file->unlink;
            }
            $self->throw("$logfile says we parsed $total_reads_parsed instead of bam meta $bam_meta->{reads} reads");
        }

        # check the fastq files are as expected
        my %extra_meta;
        foreach my $out_file (@out_files) {

            # get read and base counts from fastqcheck
            my ($these_reads, $these_bases) = (0, 0);
            my $path = $out_file->path;
            my $pipe = "$fastqcheck_exe $path |";
            open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
            while (<$fh>) { 
                # should be first line
                unless (/(\S+) sequences, (\S+) total length/) {
                    foreach my $out_file (@out_files) {
                        $out_file->unlink;
                    }
                    $self->throw("Could not get reads and bases from $pipe : $_");
                }
                ($these_reads, $these_bases) = ($1,$2);
                last;
            }
            close($fh);
            
            my $fq_meta = $out_file->metadata;
            # for expected_reads we do not use the 'reads' metadata since that
            # is based on bamcheck forward/reverse reads, which is based on
            # flags, so if flags are wrong, or reads have been removed from
            # input bam, either way resulting in mismatch of forward and reverse
            # reads, we will have discarded some reads.

            $out_file->update_stats_from_disc;
            my $output_lines = $out_file->lines;
            unless ($output_lines == $these_reads) {
                foreach my $out_file (@out_files) {
                    $out_file->unlink;
                }
                $self->throw("$path record count $output_lines disagrees with fastqcheck count $these_reads");
            }
            $extra_meta{ $out_file->id }->{reads} = $these_reads;
            $extra_meta{ $out_file->id }->{bases} = $these_bases;
        }
        
        # add/correct metadata
        foreach my $out_file (@out_files) {
            my $extra = $extra_meta{ $out_file->id } || next;

            $out_file->add_metadata({
                    bases           => $extra->{bases},
                    reads           => $extra->{reads},
                    avg_read_length => sprintf("%0.2f", $extra->{bases} / $extra->{reads})
                },
                replace_data => 1
            );
        }

        sub add_or_rm_fastq {
            # Create additional files if necessary, or get rid of the unneeded empty files
            my ($self,$path,$based_on,$out_files) = @_;

            my $pipe = "wc -l $path |";
            my $file_size;
            open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
            while (<$fh>) { 
                unless (/(\S+) $path/) {
                    foreach my $out_file (@$out_files) {
                        $out_file->unlink;
                    }
                    $self->throw("Couldn't get count from $pipe : $_");
                }
                $file_size = $1;
                last;
            }
            if ($file_size == 0 ) {
                unless (system("rm $path") == 0) {
                    foreach my $out_file (@$out_files) {
                        $out_file->unlink;
                    }
                    $self->throw("failed to rm $path");
                }
                return;
            }

            # figure out what stepstate we are for and add the fastq as an output file
            my @existing_stepoutputfiles = VRPipe::StepOutputFile->search({ file => VRPipe::File->get(path => $based_on)->id, output_key => 'fastq_files'});
            my %stepstates = map { $_->id => $_ } @existing_stepoutputfiles;
            if (keys %stepstates != 1) {
                foreach my $out_file (@$out_files) {
                    $out_file->unlink;
                }
                $self->throw("Could not get unique stepstate for $based_on");
            }

            my $fastq_meta = VRPipe::File->get(path => $based_on)->metadata;
            if ( $path =~ /_M.fastq/) {
                $fastq_meta->{paired} = 0;
                delete $fastq_meta->{mate};
            }
            elsif ($path =~ /_1.fastq/) {
                $fastq_meta->{paired} = 1;
                my $mate = $path;
                $mate =~ s/_1.fastq/_2.fastq/;
                $fastq_meta->{mate} = $mate;
            }
            else {
                $fastq_meta->{paired} = 2;
                my $mate = $path;
                $mate =~ s/_2.fastq/_1.fastq/;
                $fastq_meta->{mate} = $mate;
            }
            my $extra_file = VRPipe::File->create(path => $path,  metadata => {%$fastq_meta});

            my $step_state = $existing_stepoutputfiles[0]->stepstate;
            VRPipe::StepOutputFile->create (
                file => $extra_file->id, 
                stepstate => $step_state, 
                output_key => 'fastq_files',
            );
            push(@$out_files, $extra_file);
        }
    }
}

1;
