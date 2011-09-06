use VRPipe::Base;

class VRPipe::Steps::bam_calculate_bq with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
                 samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools'),
                 calmd_options => VRPipe::StepOption->get(description => 'options to samtools calmd',
                                                          optional => 1,
                                                          default_value => '-Erb') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                            max_files => -1, 
                                                            description => '1 or more bam files',
                                                            metadata => {reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired reads were mapped; 1=paired reads were mapped'}) };
        
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $samtools = $options->{samtools_exe};
            my $calmd_opts = $options->{calmd_options};
            if ($calmd_opts =~ /calmd/) {
                $self->throw("calmd_options should not include the calmd subcommand");
            }
            my $cmd = $bwa_exe.' calmd '.$calmd_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'samtools', version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'), summary => "samtools calmd $calmd_opts \$bam_file \$ref > \$bq_bam_file"));
            
            foreach my $bam (@{$self->input->{bam_files}}) {
                my $bam_base = $bam->basename;
                my $bq_base = $bam_base;
                $bq_base =~ s/bam$/calmd.bam/;
                my $bam_meta = $bam->metadata;
                my $bq_bam_file = $self->output_file(output_key => 'bq_bam_files', basename => $bq_base, type => 'bam', metadata => $bam_meta);
                my $bam_path = $bam->path;
                my $bq_bam_path = $bq_bam_file->path;
                my $this_cmd = "$samtools calmd $calmd_opts $bam_path $ref > $bq_bam_path";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_calculate_bq', 'calmd_and_check', [$this_cmd, $req, {output_files => [$bq_bam_file]}]);
            }
        };
    }
    method outputs_definition {
        return { bq_bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                               max_files => -1, 
                                                               description => 'a bam file with BQ tag and good NM & MD tags',
                                                               metadata => {reads => 'total number of reads (sequences)',
                                                                            paired => '0=unpaired reads were mapped; 1=paired reads were mapped'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Corrects NM & MD tags and calculates BQ (BAQ scores) to aid downstream variant calling";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method calmd_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_bam_path, $bq_bam_path) = $cmd_line =~ /(\S+) \S+ > (\S+)/;
        $in_bam_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $bq_bam_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_bam_file = VRPipe::File->get(path => $in_bam_path);
        my $bq_bam_path = VRPipe::File->get(path => $bq_bam_path);
        
        my $expected_reads = $in_bam_file->metadata->{reads};
        
        $in_bam_file->disconnect;
        $bq_bam_path->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $bq_bam_path->update_stats_from_disc(retries => 3);
        my $actual_reads = $bam_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $realigned_bam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the bq bam file, yet there were $expected_reads reads in the original bam file");
        }
    }
}

1;
