use VRPipe::Base;

class VRPipe::Steps::mpileup_bcf with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to samtools executable', optional => 1, default_value => 'samtools'),
                 samtools_mpileup_options => VRPipe::StepOption->get(description => 'samtools mpileup options excluding -f', optional => 1, default_value => '-DSV -C50 -m2 -F0.0005 -d 10000 -g'),
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to reference genome fasta'),
        };
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
            my $reference_fasta = $options->{reference_fasta};
            
            my $bam_list;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                $bam_list .= "$bam_path ";
            }
            my $bcf_file = $self->output_file(output_key => 'bcf_files', basename => "mpileup.bcf", type => 'bin');
            my $bcf_path = $bcf_file->path;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $cmd = qq[$samtools mpileup $mpileup_opts -f $reference_fasta $bam_list > $bcf_path];
            $self->dispatch_wrapped_cmd('VRPipe::Steps::mpileup_bcf', 'mpileup_bcf_and_check', [$cmd, $req, {output_files => [$bcf_file]}]); 
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'samtools', 
                                   version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'), 
                                   summary => "samtools mpileup $mpileup_opts -f \$reference_fasta \$bam_files > \$bcf_file"));
        };
    }
    method outputs_definition {
        return { bcf_files => VRPipe::StepIODefinition->get(type => 'bin', max_files => -1, description => 'a .bcf file for each set of input bam files') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run samtools mpileup for one or more bams, generating one bcf file per set of bams";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method mpileup_bcf_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($bcf_path) = $cmd_line =~ /> (\S+)$/;
        $bcf_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $bcf_file = VRPipe::File->get(path => $bcf_path);
        
        $bcf_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $bcf_file->update_stats_from_disc(retries => 3);
        
        my $correct_magic = [qw(037 213 010 004 000 000 000 000 000 377 006 000 102 103 002)];
        
        if ($bcf_file->check_magic($bcf_file->path, $correct_magic)) {
            return 1;
        }
        else {
            $bcf_file->unlink;
            $self->throw("cmd [$cmd_line] failed because index file $bcf_path had the incorrect magic");
        }
    }
}

1;
