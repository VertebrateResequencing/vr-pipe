use VRPipe::Base;

class VRPipe::Steps::bam_name_sort with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                            max_files => -1, 
                                                            description => '1 or more bam files') };
        
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            
            my $samtools = $options->{samtools_exe};
            my $req = $self->new_requirements(memory => 3000, time => 2);
            my $memory = $req->memory;
            
            # my $m = (($memory * 1000000) / 100) * 95; # we no longer use -m because in some cases samtools can use way more than the figure we specify, so safest to go with the small default
            my $opts = "sort -n";
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $in_base = $bam->basename;
                my $out_base = $in_base;
                $out_base =~ s/\.bam$//;
                $out_base .= '.name_sorted.bam';
                my $bam_meta = $bam->metadata;
                my $sort_bam_file = $self->output_file(output_key => 'name_sorted_bam_files', basename => $out_base, type => 'bam', metadata => $bam_meta);
                
                my $out_prefix = $sort_bam_file->path;
                $out_prefix =~ s/\.bam$//;
                my $this_cmd = "$samtools $opts ".$bam->path." $out_prefix";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_name_sort', 'sort_and_check', [$this_cmd, $req, {output_files => [$sort_bam_file]}]);
            }
        };
    }
    method outputs_definition {
        return { name_sorted_bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                                        max_files => -1, 
                                                                        description => 'a name-sorted bam file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Name-sorts a bam file, producing a new bam file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method sort_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /(\S+) (\S+)$/;
        $in_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path .= '.bam';
        
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
