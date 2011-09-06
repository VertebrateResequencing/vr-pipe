use VRPipe::Base;

class VRPipe::Steps::bam_realignment_around_known_indels with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
                 gatk_path => VRPipe::StepOption->get(description => 'path to GATK jar files'
                                                      optional => 1,
                                                      default_value => "$ENV{GATK}"),
                 known_indels => VRPipe::StepOption->get(description => 'ref to an array of absolute paths to VCF files of known indel sites'),
                 dbsnp => VRPipe::StepOption->get(description => 'absolute path to dbSNP VCF file'),
                 gatk_realign_options => VRPipe::StepOption->get(description => 'command line options for GATK IndelRealigner', optional => 1, default_value => '-LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing'),
                 tmp_dir => VRPipe::StepOption->get(description => 'location for tmp directories', optional => 1, default_value => ''),
                };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                            max_files => -1, 
                                                            description => '1 or more bam files', 
                                                            metadata => {reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired reads were mapped; 1=paired reads were mapped'}),
                 intervals_file => VRPipe::StepIODefinition->get(type => 'txt', description => 'intervals file created by GATK RealignerTargetCreator')};
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $gatk_exe = File::Spec->catfile($options{gatk_path}, 'GenomeAnalysisTK.jar');
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $known_indels = '';
            foreach my $indel (@{$options->{known_indels}}) {
                my $indel_file = Path::Class::File->new($indel);
                $self->throw("known indels must be given as absolute paths") unless $indel_file->is_absolute;
                $known_indels .= " -B:indels,VCF $indel_file";
            }
            
            my $dbsnp = Path::Class::File->new($options->{dbsnp});
            $self->throw("dbsnp must be an absolute path") unless $dbsnp->is_absolute;
            $dbsnp = "-B:dbsnp,VCF $dbsnp";
            
            my $intervals = $self->inputs->{intervals_file};
            
            my $realign_opts = $options->{gatk_realign_options};
            if ($realign_opts =~ /$ref|$dbsnp|IndelRealigner/) {
                $self->throw("gatk_realign_options should not include the reference, dbsnp file, known indel files or IndelRealigner task command");
            }
            
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
        return { realigned_bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                                      max_files => -1, 
                                                                      description => 'a name-sorted uncompressed bam file with improved alignments near indels',
                                                                      metadata => {reads => 'total number of reads (sequences)',
                                                                                   paired => '0=unpaired reads were mapped; 1=paired reads were mapped'}) };
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
        my $actual_reads = $bam_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $realigned_bam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the realigned bam file, yet there were $expected_reads reads in the original bam file");
        }
    }
}

1;
