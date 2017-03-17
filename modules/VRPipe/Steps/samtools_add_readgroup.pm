
=head1 NAME

VRPipe::Steps::samtools_add_readgroup - a step

=head1 DESCRIPTION

Run samtools addreplacerg to add or replace readgroup information  using
information in metadata. Currently only input files with one  @RG line are
supported. Input files must have the associated lane  metadata that will be
used as readgroup identifier, i.e. @RG ID:lane

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

class VRPipe::Steps::samtools_add_readgroup with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return {
            samtools_exe                   => VRPipe::StepOption->create(description => 'path to samtools executable',                                                                                                            optional => 1, default_value => 'samtools'),
            samtools_addreplacerg_options  => VRPipe::StepOption->create(description => 'options to samtools addreplacerg excluding -r and -R; set to "-O cram" to output CRAM',                                                  optional => 1, default_value => '-O bam'),
            readgroup_sm_from_metadata_key => VRPipe::StepOption->create(description => 'The SM of the readgroup will come from metadata associated with the CRAM; this option chooses which metadata key to get the value from', optional => 1, default_value => 'sample'),
            readgroup_ds_from_metadata_key => VRPipe::StepOption->create(description => 'The DS of the readgroup will come from metadata associated with the CRAM; this option chooses which metadata key to get the value from', optional => 1, default_value => 'study'),
            readgroup_lb_from_metadata_key => VRPipe::StepOption->create(description => 'The LB of the readgroup will come from metadata associated with the CRAM; this option chooses which metadata key to get the value from', optional => 1, default_value => 'library')
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'BAM or CRAM files',
                metadata    => {
                    lane          => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    library       => 'library name',
                    sample        => 'sample name',
                    center_name   => 'center name',
                    platform      => 'sequencing platform, eg. ILLUMINA|LS454|SOLID',
                    study         => 'name of the study',
                    platform_unit => 'platform sequencing unit',
                    reads         => 'total number of reads (sequences)',
                    optional      => ['library', 'platform_unit', 'sample', 'center_name', 'platform', 'study']
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
            my $self        = shift;
            my $options     = $self->options;
            my $samtools    = $options->{samtools_exe};
            my $arrg_opts   = $options->{samtools_addreplacerg_options};
            my $sample_key  = $options->{readgroup_sm_from_metadata_key};
            my $study_key   = $options->{readgroup_ds_from_metadata_key};
            my $library_key = $options->{readgroup_lb_from_metadata_key};
            
            my $suffix     = $arrg_opts =~ m/-O\s*cram/ ? 'cram'      : 'bam';
            my $suffix_idx = $arrg_opts =~ m/-O\s*cram/ ? 'cram.crai' : 'bam.bai';
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => qq[samtools addreplacerg -r \$rgline $arrg_opts \$input_file \$output_file]
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $file (@{ $self->inputs->{bam_files} }) {
                my $basename = $file->basename;
                $basename =~ s/\.(cr|b)am$//;
                my $meta              = $file->metadata;
                my $rg_added_bam_file = $self->output_file(
                    output_key => 'rg_added_bam_files',
                    basename   => qq[$basename.$suffix],
                    type       => "aln",
                    metadata   => $meta
                );
                my $rg_added_bam_index_file = $self->output_file(
                    output_key => 'rg_added_index_files',
                    basename   => qq[$basename.$suffix_idx],
                    type       => "bin"
                );
                $self->output_file(basename => qq[$basename.header], type => 'txt', temporary => 1);
                my $rg_line    = $self->rg_line_from_metadata($file, sample_key => $sample_key, study_key => $study_key, library_key => $library_key);
                my $input_bam  = $file->path;
                my $output_bam = $rg_added_bam_file->path;
                my $this_cmd   = "use VRPipe::Steps::samtools_add_readgroup; VRPipe::Steps::samtools_add_readgroup->add_rg_and_check(in_path => q[$input_bam], out_path => q[$output_bam], samtools => q[$samtools], arrg_opts => q[$arrg_opts], rg_line => q[$rg_line]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$rg_added_bam_file, $rg_added_bam_index_file] });
            }
        };
    }
    
    method outputs_definition {
        return {
            rg_added_bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'BAM or CRAM files with readgroup info added or replaced',
                metadata    => { reads => 'total number of reads (sequences)' }
            ),
            rg_added_index_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => -1,
                description => 'output BAI or CRAI index files',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Adds readgroup information to CRAM files using samtools addreplacerg";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method rg_line_from_metadata (ClassName|Object $self: VRPipe::File $file!, Str :$sample_key = 'sample', Str :$study_key = 'study', Str :$library_key = 'library') {
        use VRPipe::Steps::samtools_add_readgroup;
        
        # @RG ID CN DS DT FO KS LB PG PI PL PM PU SM
        my $meta = $file->metadata;
        # if (exists $tech_to_platform{ $meta->{platform} }) {
        #     $meta->{platform} = $tech_to_platform{ $meta->{platform} };
        # }
        my %rg_info = ();
        if ($file->s) {
            my $pars = VRPipe::Parser->create($file->type, { file => $file });
            %rg_info = $pars->readgroup_info();
            $pars->close;
        }
        
        my $lane    = $self->command_line_safe_string($meta->{lane});
        my $rgline  = "\@RG\\tID:" . $meta->{lane};
        my $library = $self->command_line_safe_string($meta->{$library_key} || $rg_info{LB} || '');
        $rgline .= "\\tLB:$library" if $library;
        my $platform = $self->command_line_safe_string($meta->{platform} || $rg_info{PL} || 'ILLUMINA');
        $rgline .= "\\tPL:$platform" if $platform;
        my $platform_unit = $self->command_line_safe_string($meta->{platform_unit} || $rg_info{PU} || $meta->{lane});
        $rgline .= "\\tPU:$platform_unit" if $platform_unit;
        my $sample = $self->command_line_safe_string($meta->{$sample_key} || $meta->{sample} || $rg_info{SM} || '');
        $rgline .= "\\tSM:$sample" if $sample;
        my $center = $self->command_line_safe_string($meta->{center_name} || $rg_info{CN} || '');
        $rgline .= "\\tCN:$center" if $center;
        my $study = $self->command_line_safe_string($meta->{$study_key} || $meta->{study} || $rg_info{DS} || '');
        $rgline .= "\\tDS:$study" if $study;
        
        foreach my $id (qw(DT FO KS PG PI PM)) {
            $rgline .= "\\t$id:" . $self->command_line_safe_string($rg_info{$id}) if $rg_info{$id};
        }
        return $rgline;
    }
    
    method add_rg_and_check (ClassName|Object $self: Str|File :$in_path!, Str|File :$out_path!, Str :$samtools!, Str :$arrg_opts!, Str :$rg_line!) {
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        # input bam should include only one @RG line
        my $rgs = `$samtools view -H $in_path | grep ^\@RG | wc -l`;
        $self->throw("$rgs \@RG lines present in the bam $in_path header. Only one \@RG line is accepted!") unless ($rgs == 1);
        
        # store the new header in temporary file
        my $dir      = $out_file->dir;
        my $basename = $in_file->basename;
        $basename =~ s/\.(cr|b)am$//;
        my $cmd = "$samtools view -H $in_path | sed 's/^\@RG.*\$/$rg_line/' > $dir/$basename.header";
        $in_file->disconnect;
        system($cmd) && $self->throw("failed to run [$cmd]");
        
        my $cmd_line = "$samtools reheader $dir/$basename.header $in_path";
        $cmd_line .= " | $samtools addreplacerg $arrg_opts -R " . $in_file->metadata->{lane} . " - " unless ($arrg_opts =~ m/-O\s*cram/);
        $cmd_line .= " > $out_path && $samtools index $out_path";
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{total_reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            $out_file->add_metadata({ reads => $actual_reads });
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output file, yet there were $expected_reads reads in the original file");
        }
    }
    
    sub command_line_safe_string {
        my ($self, $str) = @_;
        # add single quotes around the string if it contains spaces
        if ($str =~ /\s/) {
            $str = qq['$str'];
        }
        return $str;
    }
}

1;
