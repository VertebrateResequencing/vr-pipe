use VRPipe::Base;

# java -Xmx4g -jar GenomeAnalysisTK.jar \
#  -l INFO \
#  -R /path/to/reference.fasta \
#  -knownSites /path/to/dbsnp135.vcf.gz \
#  -I<realigned.bam>   \
#  -T CountCovariates \
#   -cov ReadGroupCovariate \
#   -cov QualityScoreCovariate \
#   -cov CycleCovariate \
#   -cov DinucCovariate \
#   -recalFile my_reads.recal_data.csv
#   -L 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y

class VRPipe::Steps::bam_count_covariates extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{$self->$orig},
                 known_sites_for_recalibration => VRPipe::StepOption->get(description => '-knownSites option(s) for GATK'),
                 gatk_count_covariates_options => VRPipe::StepOption->get(description => 'command line options for GATK CountCovariates -- must include -cov options; excludes the -knownSites option(s) which are set by another StepOption', optional => 1, default_value => '-l INFO -L 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate'),
                };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files')};
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $known_sites = $options->{known_sites_for_recalibration};
            $self->throw('No knownSites supplied') unless ($known_sites =~ /-knownSites \S+/);
            
            my $covariates_options = $options->{gatk_count_covariates_options};
            if ($covariates_options =~ /$ref|knownSites|CountCovariates/) {
                $self->throw("gatk_count_covariates_options should not include the reference, knownSites option or CountCovariates task command");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'GenomeAnalysisTK', 
                                   version => $self->gatk_version(),
                                   summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T CountCovariates -R $reference_fasta -I $bam_file -recalFile $bam_file.recal_data.csv -knownSites $known_sites_file(s) '.$covariates_options));
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_base = $bam->basename;
                my $bam_meta = $bam->metadata;
                my $recal_base = $bam_base;
                $recal_base =~ s/bam$/recal_data.csv/;
                my $recal_file = $self->output_file(output_key => 'bam_recalibration_files',
                                                  basename => $recal_base,
                                                  type => 'txt',
                                                  metadata => { source_bam => $bam->path->stringify });
                
                my $temp_dir = $options->{tmp_dir} || $recal_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe.qq[ $jvm_args -jar ].$self->jar.qq[ -T CountCovariates -R $ref -I ].$bam->path.qq[ -recalFile ].$recal_file->path.qq[ $known_sites $covariates_options];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_count_covariates', 'count_covariates_and_check', [$this_cmd, $req, {output_files => [$recal_file]}]); 
            }
        };
    }
    method outputs_definition {
        return { bam_recalibration_files => VRPipe::StepIODefinition->get(type => 'txt', 
                                                                      max_files => -1, 
                                                                      description => 'recalibration file from CountCovariates',
                                                                      metadata => { source_bam => 'the bam file used to generate this recal file' } ) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Counts covariates in a bam file creating a file to be used to recalibrate quality scores";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method count_covariates_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($recal_path) = $cmd_line =~ /-recalFile (\S+)/;
        $recal_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $recal_file = VRPipe::File->get(path => $recal_path);
        
        $recal_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $recal_file->update_stats_from_disc(retries => 3);
        my $eof = `tail -1 $recal_path`;
        chomp $eof;
        
        if ($eof eq "EOF") {
            return 1;
        }
        else {
            $recal_file->unlink;
            $self->throw("cmd [$cmd_line] there was no EOF line at the end of recalibration file $recal_path");
        }
    }
}

1;
