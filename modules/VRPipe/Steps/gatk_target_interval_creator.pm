use VRPipe::Base;

class VRPipe::Steps::bam_realignment_around_known_indels with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
                 gatk_path => VRPipe::StepOption->get(description => 'path to GATK jar files'
                                                      optional => 1,
                                                      default_value => "$ENV{GATK}/GenomeAnalysisTK.jar"),
                 known_indels => VRPipe::StepOption->get(description => 'ref to an array of absolute paths to VCF files of known indel sites'),
                 dbsnp => VRPipe::StepOption->get(description => 'absolute path to dbSNP VCF file'),
                 tmp_dir => VRPipe::StepOption->get(description => 'location for tmp directories', optional => 1, default_value => ''),
                };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $gatk_exe = $options{gatk_path};
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $bindings = $options->{realign_bindings};
            # check bindings okay
            
            my $intervals = $self->inputs->{intervals_file};
            
            my $realign_opts = $options->{gatk_realign_options};
            if ($realign_opts =~ /$ref|-o|-I|IndelRealigner/) {
                $self->throw("gatk_realign_options should not include the reference");
            }
            
            my $req = $self->new_requirements(memory => 3999, time => 2);
            my $memory = $req->memory;
            my $java_mem = int($memory * 0.9);
            my $xss = 280;
            if ($java_mem > 1000) {
                $xss = "-Xss${xss}m";
            }
            else {
                $xss = ''; # login node with small memory limit doesn't like Xss option at all
            }
            my $temp_dir = $options->{tmp_dir};
            $temp_dir = $self->tempdir($temp_dir ? (DIR => $temp_dir) : ());
            
            my $task_exe = "java -Xmx${java_mem}m -Xms${java_mem}m $xss -Djava.io.tmpdir=$temp_dir -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -jar $gatk_exe -T RealignerTargetCreator";
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'GenomeAnalysisTK', 
                                   version => VRPipe::StepCmdSummary->determine_version(qq[java -jar $gatk_exe -h], 'v([\d\.\-]+)'), 
                                   summary => 'java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -o $intervals_file '));
            
            my $intervals_file = $self->output_file(output_key => 'intervals_file',
                                              # basename => $bam->basename,
                                              type => 'txt',
                                              metadata => {source_bam => $bam->path->stringify,
                                                           reference_fasta => $ref->stringify,
                                                           reads => $bam->metadata->{reads}});
            my $this_cmd = qq[$task_exe -R $ref $bindings -o $intervals_file->path $realign_opts];
            $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_realignment_around_known_indels', 'realign_and_check', [$this_cmd, $req, {output_files => [$realigned_bam_file]}]); 
        };
    }
    method outputs_definition {
        return { intervals_file => VRPipe::StepIODefinition->get(type => 'txt', description => 'intervals file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Realigns reads around known indels to improve subsequent variant calling, producing a name-sorted uncompressed bam";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method realign_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_bam_path, $realigned_bam_path) = $cmd_line =~ /-I (\S+) -o (\S+)/;
        $in_bam_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $realigned_bam_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_bam_file = VRPipe::File->get(path => $in_bam_path);
        my $realigned_bam_file = VRPipe::File->get(path => $realigned_bam_path);
        
        my $expected_reads = $in_bam_file->metadata->{reads};
        
        $in_bam_file->disconnect;
        $realigned_bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $realigned_bam_file->update_stats_from_disc(retries => 3);
        my $actual_reads = $bam_file->lines;
        
        if ($actual_lines == $expected_lines) {
            return 1;
        }
        else {
            $realigned_bam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_lines lines were generated in the realigned bam file, yet there were $expected_lines lines in the original bam file");
        }
    }
}

1;
