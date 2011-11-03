use VRPipe::Base;

class VRPipe::Steps::bam_add_readgroup extends VRPipe::Steps::picard {
    around options_definition {
        return { %{$self->$orig},
                 picard_add_readgroups_options => VRPipe::StepOption->get(description => 'options for picard AddOrReplaceReadGroups', optional => 1, default_value => 'VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0'),
                };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                            max_files => -1,
                                                            description => 'bam files to have header replaced',
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                         library => 'library name',
                                                                         sample => 'sample name',
                                                                         center_name => 'center name',
                                                                         platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                         study => 'name of the study',
                                                                         platform_unit => 'platform sequencing unit',,
                                                                         reads => 'total number of reads (sequences)',
                                                                         optional => ['library', 'platform_unit', 'sample', 'center_name', 'platform', 'study']}) };
    }
    method body_sub {
        return sub {
            use VRPipe::Utils::picard;
            
            my $self = shift;
            my $options = $self->options;
            my $picard = VRPipe::Utils::picard->new(picard_path => $options->{picard_path}, java_exe => $options->{java_exe});
            my $picard_jar = Path::Class::File->new($picard->picard_path, 'AddOrReplaceReadGroups.jar');
            
            my $opts = $options->{picard_add_readgroups_options};
            my $rginfo = 'RGID=$lane RGLB=$library RGPL=$platform RGPU=$platform_unit RGSM=$sample RGCN=$centre RGDS=$study';
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'picard', 
                                   version => $picard->determine_picard_version(),
                                   summary => 'java $jvm_args -jar AddOrReplaceReadGroups.jar INPUT=$bam_file OUTPUT=$rg_added_bam_file'." $rginfo $opts"));
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $memory = $req->memory;
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $meta = $bam->metadata;
                my $rginfo_cmd;
                $rginfo_cmd = "RGID=".$meta->{lane};
                my $library = $meta->{library} || 'unknown_library';
                $rginfo_cmd .= " RGLB=$library";
                my $platform = $meta->{platform} || 'unknown_platform';
                $rginfo_cmd .= " RGPL=$platform";
                my $platform_unit = $meta->{platform_unit} || 'unknown_platform_unit';
                $rginfo_cmd .= " RGPU=$platform_unit";
                my $sample = $meta->{sample} || 'unknown_sample';
                $rginfo_cmd .= " RGSM=$sample";
                my $centre = $meta->{centre_name} || 'unknown_centre';
                $rginfo_cmd .= " RGCN=$centre";
                my $study = $meta->{study} || 'unknown_study';
                $rginfo_cmd .= " RGDS=$study";
                
                my $rg_added_bam_file = $self->output_file(output_key => 'rg_added_bam_files',
                                                  basename => $bam->basename,
                                                  type => 'bam',
                                                  metadata => $meta);
                
                my $temp_dir = $options->{tmp_dir} || $rg_added_bam_file->dir;
                my $jvm_args = $picard->jvm_args($memory, $temp_dir);
                
                my $this_cmd = $picard->java_exe." $jvm_args -jar $picard_jar INPUT=".$bam->path." OUTPUT=".$rg_added_bam_file->path." $rginfo_cmd $opts";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_add_readgroup', 'add_rg_and_check', [$this_cmd, $req, {output_files => [$rg_added_bam_file]}]); 
            }
        };
    }
    method outputs_definition {
        return { rg_added_bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                                  max_files => -1,
                                                                  description => 'uncompressed bam files with readgroup info added',
                                                                  metadata => {reads => 'total number of reads (sequences)'}) };
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
        $in_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
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
}

1;
