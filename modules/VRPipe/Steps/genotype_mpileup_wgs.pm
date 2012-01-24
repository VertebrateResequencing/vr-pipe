use VRPipe::Base;

class VRPipe::Steps::genotype_mpileup_wgs with VRPipe::StepRole {

	method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
        		 samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools'),
                 samtools_options => VRPipe::StepOption->get(description => 'options provided for samtools mpileup', 
                                                             optional => 1, 
                                                             default_value => '-ugDI -d 1000')
        };
    }

    method inputs_definition {
        return { 
        		 bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'bam files for bcf production'),
        		 hapmap_file => VRPipe::StepIODefinition->get(type => 'txt', description => 'text file with the hapmap output from the snp binary file'),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            my $samtools = $options->{samtools_exe};
            my $samtools_options = $options->{samtools_options};
            my $hapmap_path = $self->inputs->{hapmap_file}->[0]->path;
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $bam_dir = $bam->dir;
                my $bam_basename = $bam->basename;
                my $sample;
                my $meta = $bam->metadata;
                if ($meta->{sample}) {
                	$sample = $meta->{sample};
                }
                else {
             	    my $bam_irods_dir = (split('_', $bam_basename))[0];
             	    my @samtools_sample = split("\t", `samtools view -H $bam_path | grep SM`);
             	    for my $sam_out ( @samtools_sample ) {
             	    	if ( $sam_out =~ s/^SM:// ) {
             	    		$sample = $sam_out;
             	    	}
             	    }
                }
                $self->throw("unable to obtain sample name for this bam file") unless $sample;
                my $bcf_file = $self->output_file(output_key => 'bcf_files_with_metadata',
                                              basename => $bam_basename.'.bcf',
                                              type => 'bin',
                                              metadata => {sample => $sample});
                my $bcf_path = $bcf_file->path;
            	my $cmd = qq[$samtools mpileup $samtools_options -l $hapmap_path -f $ref $bam_path > $bcf_path];
                $self->dispatch([$cmd, $req, {output_files => [$bcf_file]}]);
            }
        return 1;
        };
    }
    
    method outputs_definition {
        return { bcf_files_with_metadata => VRPipe::StepIODefinition->get(type => 'bin',
                                                                  max_files => -1,
                                                                  description => 'bcf file produced for bam file using samtools',
                                                                  metadata => {sample => 'sample name'}
                                                                  )
         };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "This step uses samtools mpileup to create a bcf file from a bam file using a hapmap file";
    }
    method max_simultaneous {
        return 0;
    }
}

1;
