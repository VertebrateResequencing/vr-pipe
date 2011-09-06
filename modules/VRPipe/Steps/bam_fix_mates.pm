use VRPipe::Base;

# net.sf.picard.sam.FixMateInformation INPUT=[/lustre/scratch106/projects/uk10k/TRACKING/UK10K_COHORT_TWINSUK/QTL211713/SLX/HUMbpgRVWDIAAPEI_11/110501_I283_FCB012RABXX_L5_HUMbpgRVWDIAAPEI-11/24175.pe.realigned.bam] OUTPUT=/lustre/scratch106/projects/uk10k/TRACKING/UK10K_COHORT_TWINSUK/
# QTL211713/SLX/HUMbpgRVWDIAAPEI_11/110501_I283_FCB012RABXX_L5_HUMbpgRVWDIAAPEI-11/24175.pe.realigned.sorted.working.bam SORT_ORDER=coordinate TMP_DIR=/lustre/scratch106/tmp/1yyuaitrS2 VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=990000    VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 CREATE_INDEX=false CREATE_M
# D5_FILE=false

class VRPipe::Steps::bam_fix_mates with VRPipe::StepRole {
    method options_definition {
        return { picard_path => VRPipe::StepOption->get(description => 'path to Picard jar files'
                                              optional => 1,
                                              default_value => "$ENV{PICARD}"),
                 tmp_dir => VRPipe::StepOption->get(description => 'location for tmp directories', optional => 1, default_value => ''),
               };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more name-sorted bam files') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $exe = File::Spec->catfile($options{picard_path}, 'FixMateInformation.jar');
            
            my $req = $self->new_requirements(memory => 4000, time => 3);
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
            
            my $task_exe = "java -Xmx${java_mem}m -Xms${java_mem}m $xss -Djava.io.tmpdir=$temp_dir -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -jar $gatk_exe -T IndelRealigner";
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'GenomeAnalysisTK', 
                                   version => VRPipe::StepCmdSummary->determine_version(qq[java -jar $gatk_exe -h], 'v([\d\.\-]+[a-z\d]*)'), 
                                   summary => 'java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -B:dbsnp,VCF $dbsnp -B:indels,VCF $known_indels -I $in_bam -o $realigned_bam -targetIntervals $intervals_file '.$realign_opts));
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_base = $bam->basename;
                my $bam_meta = $bam->metadata;
                my $realigned_base = $bam_base;
                $realigned_base =~ s/bam$/realign.bam/;
                my $realigned_bam_file = $self->output_file(output_key => 'realigned_bam_files',
                                                  basename => $realigned_base,
                                                  type => 'bam',
                                                  metadata => $bam_meta);
                my $this_cmd = qq[$task_exe -R $ref $dbsnp $known_indels -I $bam->path -o $realigned_bam_file->path -targetIntervals $intervals->path $realign_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_realignment_around_known_indels', 'realign_and_check', [$this_cmd, $req, {output_files => [$realigned_bam_file]}]); 
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
}

1;
