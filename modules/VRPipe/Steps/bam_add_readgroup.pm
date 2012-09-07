
=head1 NAME

VRPipe::Steps::bam_add_readgroup - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::bam_add_readgroup extends VRPipe::Steps::picard {
    use VRPipe::Parser;
    
    around options_definition {
        return {
            %{ $self->$orig },
            picard_add_readgroups_options  => VRPipe::StepOption->create(description => 'options for picard AddOrReplaceReadGroups',                                                                                             optional => 1, default_value => 'VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0'),
            readgroup_sm_from_metadata_key => VRPipe::StepOption->create(description => 'The SM of the readgroup will come from metadata associated with the bam; this option chooses which metadata key to get the value from', optional => 1, default_value => 'sample')
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'bam files to have header replaced',
                metadata    => {
                    lane          => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    library       => 'library name',
                    sample        => 'sample name',
                    center_name   => 'center name',
                    platform      => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study         => 'name of the study',
                    platform_unit => 'platform sequencing unit',
                    reads         => 'total number of reads (sequences)',
                    optional      => ['library', 'platform_unit', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            my $picard_jar = $self->jar('AddOrReplaceReadGroups.jar');
            
            my $sample_key = $options->{readgroup_sm_from_metadata_key};
            my $opts       = $options->{picard_add_readgroups_options};
            my $rginfo     = 'RGID=$lane RGLB=$library RGPL=$platform RGPU=$platform_unit RGSM=$' . $sample_key . ' RGCN=$centre RGDS=$study';
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'picard',
                    version => $self->picard_version(),
                    summary => 'java $jvm_args -jar AddOrReplaceReadGroups.jar INPUT=$bam_file OUTPUT=$rg_added_bam_file' . " $rginfo $opts"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $memory = $req->memory;
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta = $bam->metadata;
                my $rginfo_cmd;
                $rginfo_cmd = "RGID=" . $meta->{lane};
                my $library = $self->command_line_safe_string($meta->{library} || 'unknown_library');
                $rginfo_cmd .= " RGLB=$library";
                my $platform = $self->command_line_safe_string($meta->{platform} || 'unknown_platform');
                $rginfo_cmd .= " RGPL=$platform";
                my $platform_unit = $self->command_line_safe_string($meta->{platform_unit} || $meta->{lane});
                $rginfo_cmd .= " RGPU=$platform_unit";
                my $sample = $self->command_line_safe_string($meta->{$sample_key} || $meta->{sample} || 'unknown_sample');
                $rginfo_cmd .= " RGSM=$sample";
                my $centre = $self->command_line_safe_string($meta->{centre_name} || 'unknown_centre');
                $rginfo_cmd .= " RGCN=$centre";
                my $study = $self->command_line_safe_string($meta->{study} || 'unknown_study');
                $rginfo_cmd .= " RGDS=$study";
                
                my $rg_added_bam_file = $self->output_file(
                    output_key => 'rg_added_bam_files',
                    basename   => $bam->basename,
                    type       => 'bam',
                    metadata   => $meta
                );
                
                my $temp_dir = $options->{tmp_dir} || $rg_added_bam_file->dir;
                my $jvm_args = $self->jvm_args($memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . " $jvm_args -jar $picard_jar INPUT=" . $bam->path . " OUTPUT=" . $rg_added_bam_file->path . " $rginfo_cmd $opts";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_add_readgroup', 'add_rg_and_check', [$this_cmd, $req, { output_files => [$rg_added_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            rg_added_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'uncompressed bam files with readgroup info added',
                metadata    => { reads => 'total number of reads (sequences)' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Adds readgroup information to bam files";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method add_rg_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /INPUT=(\S+) OUTPUT=(\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        # first check that it is necessary to do anything
        my $pars = VRPipe::Parser->create('bam', { file => $in_path });
        my %readgroup_info = $pars->readgroup_info();
        $in_file->disconnect;
        if (keys %readgroup_info == 1) {
            $pars->get_fields('RG');
            $pars->next_record();
            $pars->close;
            my $record_rg = $pars->parsed_record()->{RG};
            if ($record_rg && defined $readgroup_info{$record_rg}) {
                my $actual_data = $readgroup_info{$record_rg};
                my $all_match   = 1;
                foreach my $tag (qw(LB PL PU SM CN DS)) {
                    my $desired_val;
                    if ($cmd_line =~ /RG$tag='([^']+)'/) {
                        $desired_val = $1;
                    }
                    else {
                        ($desired_val) = $cmd_line =~ / RG$tag=(\S+)/;
                    }
                    next if $desired_val eq $actual_data->{$tag};
                    $all_match = 0;
                    last;
                }
                
                if ($all_match) {
                    $in_file->move($out_file);
                    return;
                }
            }
        }
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output bam file, yet there were $expected_reads reads in the original bam file");
        }
    }
    
    method command_line_safe_string (Str $str) {
        # add single quotes around the string if it contains spaces
        if ($str =~ /\s/) {
            $str = qq['$str'];
        }
        return $str;
    }
}

1;
