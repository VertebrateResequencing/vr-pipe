use VRPipe::Base;

class VRPipe::Steps::mpileup_bcf_hapmap extends VRPipe::Steps::mpileup_bcf {
#for wgs genotyping override samtools_mpileup_options with '-ugDI -d 1000' 
   
   method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'bam files for bcf production'),
        		 hapmap_file => VRPipe::StepIODefinition->get(type => 'txt', description => 'hapmap sites file for -l option'),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            use VRPipe::Steps::mpileup_bcf;
            
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            my $samtools = $options->{samtools_exe};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            my $hapmap_path = $self->inputs->{hapmap_file}->[0]->path;
            $mpileup_opts .= " -l $hapmap_path";
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
            	my $cmd = qq[$samtools mpileup $mpileup_opts -f $ref $bam_path > $bcf_path];
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
    
    method description {
        return "Run samtools mpileup for each bam using the hapmap sites file provided to generate one bcf file with associated sample name metadata per bam";
    }
}

1;
