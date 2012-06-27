use VRPipe::Base;

class VRPipe::Steps::mpileup_bcf with VRPipe::StepRole {
    method options_definition {
	return { samtools_exe => VRPipe::StepOption->get(description => 'path to samtools executable', 
							 optional => 1, 
							 default_value => 'samtools'),
		 samtools_mpileup_options => VRPipe::StepOption->get(description => 'samtools mpileup options excluding -f', 
								     optional => 1, 
								     default_value => '-C50 -aug'),
		 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to reference genome fasta'),
		 max_cmdline_bams => VRPipe::StepOption->get(description => 'max number of bam filenames to allow on command line',
							     optional => 1,
							     default_value => 10),
		 interval_list => VRPipe::StepOption->get(description => 'absolute path to targets interval list file for -l option', 
							  optional => 1),
		 snp_coordinates_file => VRPipe::StepOption->get(description => 'absolute path to file containing the coordinates of snp locations',
		                                                 optional => 1),					
		 mimimum_calls => VRPipe::StepOption->get(description => 'minumum expected number of variant calls',
							  optional => 1,
							  default_value => 0),
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
	    
            my $interval_list = $options->{interval_list};
	    $mpileup_opts .= " -l $interval_list " if $interval_list;
	    
	    my $max_cmdline_bams = $options->{max_cmdline_bams};
	    if (scalar (@{$self->inputs->{bam_files}}) > $max_cmdline_bams) {
		$self->warn("[todo] Generate a bam fofn");
	    }
	    
	    my $bam_list;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
		$bam_list .= "$bam_path ";
            }
	    my $bcf_file = $self->output_file(output_key => 'bcf_files', basename => "mpileup.bcf", type => 'bin');
	    my $bcf_path = $bcf_file->path;
	    
            my $req = $self->new_requirements(memory => 500, time => 1);
	    my $cmd = qq[$samtools mpileup $mpileup_opts -f $reference_fasta $bam_list > $bcf_path];
	    $self->dispatch([$cmd, $req, {output_files => [$bcf_file]}]); 
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
}

1;
