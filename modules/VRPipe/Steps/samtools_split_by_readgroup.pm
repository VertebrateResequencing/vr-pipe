
=head1 NAME

VRPipe::Steps::samtools_split_by_readgroup - a step

=head1 DESCRIPTION

Run samtools split to split a BAM or CRAM file into BAM or CRAM files
per-readgroup

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

class VRPipe::Steps::samtools_split_by_readgroup extends VRPipe::Steps::samtools_add_readgroup {
    use VRPipe::Parser;
    
    method options_definition {
        return {
            samtools_exe           => VRPipe::StepOption->create(description => 'path to samtools executable',                                          optional => 1, default_value => 'samtools'),
            samtools_split_options => VRPipe::StepOption->create(description => 'options to samtools split excluding; set to "-O cram" to output CRAM', optional => 1, default_value => ''),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'BAM or CRAM files to be split by readgroup',
                metadata    => {
                    lane          => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    library       => 'library name',
                    sample        => 'sample name',
                    center_name   => 'center name',
                    platform      => 'sequencing platform, eg. ILLUMINA|LS454|SOLID',
                    study         => 'name of the study',
                    platform_unit => 'platform sequencing unit',
                    optional      => ['lane', 'library', 'platform_unit', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    # CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and PACBIO
    our %tech_to_platform = (
        SLX       => 'ILLUMINA',
        '454'     => 'LS454',
        SOLID     => 'SOLID',
        ILLUMINA  => 'ILLUMINA',
        'LS454'   => 'LS454',
        ABI_SOLID => 'SOLID',
    );
    
    method body_sub {
        return sub {
            my $self       = shift;
            my $options    = $self->options;
            my $samtools   = $options->{samtools_exe};
            my $split_opts = $options->{samtools_split_options};
            
            my $suffix = $split_opts =~ m/-O\s*cram/ ? 'cram' : 'bam';
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => 'samtools split -u $orphan_file:$orphan_hdr -f %*_%!.bam $input_file'
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $memory = $req->memory;
            
            my $idx = 1;
            foreach my $aln (@{ $self->inputs->{bam_files} }) {
                my $meta           = $aln->metadata;
                my $parser         = VRPipe::Parser->create($aln->type, { file => $aln });
                my %readgroup_info = $parser->readgroup_info();
                my $basename       = $aln->basename;
                $basename =~ s/\.(cr|b)am$//;
                my @outfiles;
                while (my ($rgid, $rg_info) = each %readgroup_info) {
                    my $rg_meta = {};
                    while (my ($key, $val) = each %$meta) {
                        $$rg_meta{qq[$key]} = $val;
                    }
                    $rg_meta->{lane} ||= $rgid;
                    if ($rg_info->{LB}) {
                        $$rg_meta{library} ||= $rg_info->{LB};
                    }
                    if ($rg_info->{SM}) {
                        $$rg_meta{sample} ||= $rg_info->{SM};
                    }
                    if ($rg_info->{PL}) {
                        $$rg_meta{platform} ||= $rg_info->{PL};
                    }
                    if ($rg_info->{PU}) {
                        $$rg_meta{platform_unit} ||= $rg_info->{PU};
                    }
                    if ($rg_info->{CN}) {
                        $$rg_meta{center_name} ||= $rg_info->{CN};
                    }
                    if ($rg_info->{DS}) {
                        $$rg_meta{study} ||= $rg_info->{DS};
                    }
                    
                    push @outfiles,
                      $self->output_file(
                        output_key => 'rg_split_files',
                        sub_dir    => $idx,
                        basename   => $basename . "_${rgid}.$suffix",
                        type       => $suffix,
                        metadata   => $rg_meta
                      );
                }
                $parser->close;
                
                # create files for unknown readgroups. will only include as
                # output_files later if file contains any reads
                my $unknown_rg_head = $self->output_file(
                    sub_dir   => $idx,
                    basename  => $basename . "_U${idx}.sam",
                    temporary => 1,
                    type      => 'txt',
                )->path;
                my $unknown_rg_path = $unknown_rg_head;
                $unknown_rg_path =~ s/sam$/$suffix/;
                my $unknown_rg_file = VRPipe::File->create(path => $unknown_rg_path);
                $unknown_rg_file->add_metadata({ %$meta, lane => $idx });
                
                my @outs = map { $_->path } @outfiles;
                my $cmd = "$samtools split -u $unknown_rg_path:$unknown_rg_head -f " . file($outfiles[0]->dir, '%*_%\!.bam')->stringify . " " . $aln->path;
                my $this_cmd = "use VRPipe::Steps::samtools_split_by_readgroup; VRPipe::Steps::samtools_split_by_readgroup->split_and_check(q[$cmd], [qw(@outs)], unknown_rg_path => q[$unknown_rg_path]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
                
                $idx++;
            }
        };
    }
    
    method outputs_definition {
        return {
            rg_split_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'BAM or CRAM files containing a single readgroup each',
                metadata    => { lane => '', reads => 'total number of reads (sequences)' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Splits BAM or CRAM files to output a BAM or CRAM files per readgroup";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method split_and_check (ClassName|Object $self: Str $cmd_line, ArrayRef[Str|File] $outpaths!, Str|File :$unknown_rg_path!) {
        my ($in_path) = $cmd_line =~ /(\S+)$/;
        $in_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file = VRPipe::File->get(path => $in_path);
        my @outfiles;
        foreach my $path (@$outpaths) {
            push @outfiles, VRPipe::File->get(path => $path);
        }
        
        # create the header file for the unknown_rg reads
        my $unknown_rg_file = VRPipe::File->get(path => $unknown_rg_path);
        my $unknown_rg_header_path = $unknown_rg_path;
        $unknown_rg_header_path =~ s/(cr|b)am$/sam/;
        my $unknown_rg_line        = $self->rg_line_from_metadata($unknown_rg_file);
        my $unknown_rg_header_file = VRPipe::File->get(path => $unknown_rg_header_path);
        my $fh                     = $unknown_rg_header_file->openw;
        foreach my $line (@{ $in_file->header_lines }) {
            next if ($line =~ /^\@RG/);
            print $fh "$line\n";
        }
        print $fh "$unknown_rg_line\n";
        $fh->close;
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $total_output_reads = 0;
        my @output_reads;
        foreach my $outfile (@outfiles) {
            $outfile->update_stats_from_disc(retries => 3);
            push @output_reads, $outfile->num_records;
            $total_output_reads += $output_reads[-1];
        }
        
        # if the unknown_rg file has records, add it as an output file
        # otherwise delete
        my $unknown_rg_reads = $unknown_rg_file->num_records;
        if ($unknown_rg_reads) {
            $total_output_reads += $unknown_rg_reads;
            push @output_reads, $unknown_rg_reads;
            push @outfiles,     $unknown_rg_file;
            my @existing_stepoutputfiles = VRPipe::StepOutputFile->search({ file => $outfiles[0]->id, output_key => 'rg_split_files' });
            my %stepstates = map { $_->id => $_ } @existing_stepoutputfiles;
            if (keys %stepstates != 1) {
                $self->throw("Could not get unique stepstate for $$outpaths[0]");
            }
            my $step_state = $existing_stepoutputfiles[0]->stepstate;
            VRPipe::StepOutputFile->create(
                file       => $unknown_rg_file->id,
                stepstate  => $step_state,
                output_key => 'rg_split_files',
            );
        }
        else {
            $unknown_rg_file->unlink;
        }
        
        if ($total_output_reads == $expected_reads) {
            foreach my $outfile (@outfiles) {
                my $reads = shift @output_reads;
                $outfile->add_metadata({ reads => $reads });
            }
            return 1;
        }
        else {
            foreach my $outfile (@outfiles) {
                $outfile->unlink;
            }
            $self->throw("cmd [$cmd_line] failed because $total_output_reads reads were generated in the output files, yet there were $expected_reads reads in the original input file");
        }
    }

}

1;
