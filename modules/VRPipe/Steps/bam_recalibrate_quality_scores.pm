use VRPipe::Base;

# java -Xmx4g -jar GenomeAnalysisTK.jar \
#  -l INFO \
#  -R /path/to/reference.fasta \
#  -I<realigned.bam>   \
#  -T TableRecalibration \
#   -o<realigned.recalibrated.bam>   \
#   -recalFile my_reads.recal_data.csv

class VRPipe::Steps::bam_recalibrate_quality_scores extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{$self->$orig},
                 bam_recalibration_options => VRPipe::StepOption->get(description => 'command line options for GATK TableRecalibration', optional => 1, default_value => '-l INFO --disable_bam_indexing'),
                };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more coordinate-sorted bam files'),
                 bam_recalibration_files => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1, description => '1 or more bam recal files from count covariates step') };
    }
    method body_sub {
        return sub {
            use VRPipe::Utils::gatk;
            
            my $self = shift;
            my $options = $self->options;
            my $gatk = VRPipe::Utils::gatk->new(gatk_path => $options->{gatk_path}, java_exe => $options->{java_exe});
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $recal_opts = $options->{bam_recalibration_options};
            if ($recal_opts =~ /$ref|recalFile|TableRecalibration/) {
                $self->throw("bam_recalibration_options should not include the reference, recalFile option or TableRecalibration task command");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'GenomeAnalysisTK', 
                                   version => $gatk->determine_gatk_version(),
                                   summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T TableRecalibration -R $reference_fasta -recalFile $bam_file.recal_data.csv -I $bam_file -o $recalibrated_bam_file '.$recal_opts));
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            foreach my $recal_file (@{$self->inputs->{bam_recalibration_files}}) {
                my $bam_path = $recal_file->metadata->{source_bam};
                my $bam = VRPipe::File->get(path => $bam_path);
                my $bam_base = $bam->basename;
                my $bam_meta = $bam->metadata;
                my $recal_base = $bam_base;
                $recal_base =~ s/bam$/recal.bam/;
                my $recal_bam_file = $self->output_file(output_key => 'recalibrated_bam_files',
                                                  basename => $recal_base,
                                                  type => 'bam',
                                                  metadata => $bam_meta);
                
                my $temp_dir = $options->{tmp_dir} || $recal_bam_file->dir;
                my $jvm_args = $gatk->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $gatk->java_exe.qq[ $jvm_args -jar ].$gatk->jar.qq[ -T TableRecalibration -R $ref -recalFile ].$recal_file->path.qq[ -I ].$bam->path.qq[ -o ].$recal_bam_file->path.qq[ $recal_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_recalibrate_quality_scores', 'recal_and_check', [$this_cmd, $req, {output_files => [$recal_bam_file]}]); 
            }
        };
    }
    method outputs_definition {
        return { recalibrated_bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                                      max_files => -1, 
                                                                      description => 'a bam file with recalibrated quality scores; OQ tag holds the original quality scores',
                                                                      ) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Recalibrate quality scores of each mapped base using GATK";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method recal_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /-I (\S+) -o (\S+)/;
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
