use VRPipe::Base;

class VRPipe::Steps::bam_to_cram extends VRPipe::Steps::cramtools {
    around options_definition {
        return { %{$self->$orig},
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to reference fasta file'),
                 cramtools_bam_to_cram_options => VRPipe::StepOption->get(description => 'Options for cramtools cram command to convert bam to cram', optional => 1),
               };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more coordinate-sorted bam files'),
                 bai_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => '1 or more bam index files') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $cramtools = VRPipe::Utils::cramtools->new(cramtools_path => $options->{cramtools_path}, java_exe => $options->{java_exe});
            my $cramtools_jar = Path::Class::File->new($cramtools->cramtools_path, 'cramtools.jar');
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $opts = $options->{cramtools_bam_to_cram_options};
            if ($opts =~ /$ref|--output-cram-file|--input-bam-file|--reference-fasta-file/) {
                $self->throw("cramtools_bam_to_cram_options should not include the reference, input or output options");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'cramtools', 
                                   version => $cramtools->determine_cramtools_version(),
                                   summary => 'java $jvm_args -jar cramtools.jar cram --input-bam-file $bam_file --output-cram-file $cram_file --reference-fasta-file $reference_fasta '.$opts));
            
            my $req = $self->new_requirements(memory => 4000, time => 3);
            my $memory = $req->memory;
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_base = $bam->basename;
                my $bam_meta = $bam->metadata;
                my $cram_base = $bam_base;
                $cram_base =~ s/bam$/cram/;
                my $cram_file = $self->output_file(output_key => 'cram_files',
                                                   basename => $cram_base,
                                                   type => 'cram',
                                                   metadata => $bam_meta);
                
                my $temp_dir = $options->{tmp_dir} || $cram_file->dir;
                my $jvm_args = $cramtools->jvm_args($memory, $temp_dir);
                
                my $this_cmd = $cramtools->java_exe." $jvm_args -jar $cramtools_jar cram --input-bam-file ".$bam->path." --output-cram-file ".$cram_file->path." --reference-fasta-file $ref $opts";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_to_cram', 'cram_and_check', [$this_cmd, $req, {output_files => [$cram_file]}]); 
            }
        };
    }
    method outputs_definition {
        return { cram_files => VRPipe::StepIODefinition->get(type => 'cram', max_files => -1, description => 'a cram file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Converts bam files to cram files using cramtools";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method cram_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /--input-bam-file (\S+) --output-cram-file (\S+)/;
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
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output cram file, yet there were $expected_reads reads in the original bam file");
        }
    }
}

1;
