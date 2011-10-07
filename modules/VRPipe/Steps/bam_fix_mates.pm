use VRPipe::Base;

#java -Djava.io.tmpdir=/path/to/tmpdir \  [this argument recommended when dealing with large input]
#   -jar FixMateInformation.jar \
#   INPUT=<input1.bam> \
#   OUTPUT=<fixedBam.bam> \
#   SO=coordinate \
#   VALIDATION_STRINGENCY=SILENT


class VRPipe::Steps::bam_fix_mates with VRPipe::StepRole {
    method options_definition {
        return { picard_fix_mates_options => VRPipe::StepOption->get(description => '', optional => 1, default_value => 'SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0'),
                 picard_path => VRPipe::StepOption->get(description => 'path to Picard jar files', optional => 1, default_value => "$ENV{PICARD}"),
                 java_exe => VRPipe::StepOption->get(description => 'path to your java executable', optional => 1, default_value => 'java'),
                 tmp_dir => VRPipe::StepOption->get(description => 'location for tmp directories; defaults to working directory', optional => 1),
               };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more name-sorted bam files') };
    }
    method body_sub {
        return sub {
            use VRPipe::Utils::picard;
            
            my $self = shift;
            my $options = $self->options;
            my $picard = VRPipe::Utils::picard->new(picard_path => $options->{picard_path}, java_exe => $options->{java_exe});
            my $fixmates_jar = Path::Class::File->new($picard->picard_path, 'FixMateInformation.jar');
            
            my $fixmate_options = $options->{picard_fix_mates_options};
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'picard', 
                                   version => $picard->determine_picard_version(),
                                   summary => 'java $jvm_args -jar FixMateInformation.jar INPUT=$bam_file OUTPUT=$fixmate_bam_file '.$fixmate_options));
            
            my $req = $self->new_requirements(memory => 4000, time => 3);
            my $memory = $req->memory;
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_base = $bam->basename;
                my $bam_meta = $bam->metadata;
                my $mate_fixed_base = $bam_base;
                $mate_fixed_base =~ s/bam$/sort.bam/;
                my $mate_fixed_file = $self->output_file(output_key => 'fixmate_bam_files',
                                                  basename => $mate_fixed_base,
                                                  type => 'bam',
                                                  metadata => $bam_meta);
                
                my $temp_dir = $options->{tmp_dir} || $mate_fixed_file->dir;
                my $jvm_args = $picard->jvm_args($memory, $temp_dir);
                
                my $this_cmd = $picard->java_exe." $jvm_args -jar $fixmates_jar INPUT=".$bam->path." OUTPUT=".$mate_fixed_file->path." $fixmate_options";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_fix_mates', 'fix_mates_and_check', [$this_cmd, $req, {output_files => [$mate_fixed_file]}]); 
            }
        };
    }
    method outputs_definition {
        return { fixmate_bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'a coordinate-sorted uncompressed bam file with fixed mates') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Fixes mate information and coordinate sorts a name-sorted bam file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method fix_mates_and_check (ClassName|Object $self: Str $cmd_line) {
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
