use VRPipe::Base;

class VRPipe::Steps::mpileup_bcf with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools'),
                 samtools_mpileup_options => VRPipe::StepOption->get(description => 'samtools mpileup options',
                                                          optional => 1,
                                                          default_value => '-C50 -aug') }; # from SNPS.pm
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files to call variants') };
    }

    method body_sub {
        return sub {

            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            if (!$mpileup_opts =~ /g/) {
                $self->throw("samtools_mpileup_options must include -g flag to generate bcf");
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{$self->inputs->{bam_files}}) {

                my $bam_path = $bam->path;
				my $basename = $bam->basename;
				$basename =~ s/bam$/bcf/;
                my $bcf_file = $self->output_file(output_key => 'bcf_files', basename => $basename, type => 'bin');
                my $bcf_path = $bcf_file->path;

                my $cmd = qq[$samtools mpileup $mpileup_opts $bam_path > $bcf_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::mpileup_bcf', 'mpileup_bcf', [$cmd, $req, {output_files => [$bcf_file]}]); 
            }
        };
    }
    method outputs_definition {
        return { bcf_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => 'a .bcf file for each input bam file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run samtools mpileup and generate one bcf file per input bam";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method mpileup_bcf (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $bcf_path) = $cmd_line =~ / (\S+) > (\S+)$/;
        
        my $bam_file = VRPipe::File->get(path => $bam_path);
        my $bcf_file = VRPipe::File->get(path => $bcf_path);
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $bcf_file->update_stats_from_disc;
        
    }

}

1;
