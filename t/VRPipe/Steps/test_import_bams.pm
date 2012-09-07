use VRPipe::Base;

class VRPipe::Steps::test_import_bams with VRPipe::StepRole {
    method options_definition {
        return {};
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => 'one or more bam files') };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $req = $self->new_requirements(memory => 50, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $local_bam = $self->output_file(
                    output_key => 'local_bam_files',
                    basename   => $bam->basename,
                    type       => 'bam',
                    metadata   => $bam->metadata
                );
                
                my ($bam_path, $local_bam_path) = ($bam->path, $local_bam->path);
                my $this_cmd = "cp $bam_path $local_bam_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::test_import_bams', 'import_and_add_metadata', [$this_cmd, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            local_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'bams imported into the hashed directories',
                metadata    => { reads => 'number of reads in the bam file' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Copy bams into hashed directories and add metadata";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method import_and_add_metadata (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /cp (\S+) (\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        $out_file->add_metadata({ reads => $actual_reads });
        
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
