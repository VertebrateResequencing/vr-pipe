use VRPipe::Base;

class VRPipe::Steps::genotype_mpileup_wgs with VRPipe::StepRole {

	method options_definition {
        return { 
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
        		 snp_binary_file => VRPipe::StepOption->get(description => 'absolute path to snp binary file used for genotyping'),
        		 samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools'),
                 bin2hapmap_exe => VRPipe::StepOption->get(description => 'path to your hapmap2bin executable',
                                                           optional => 1,
                                                           default_value => 'bin2hapmap')
        };
    }

    method inputs_definition {
        return { 
        		 bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => 'bam files for bcf production')
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            my $snp_file = Path::Class::File->new($options->{snp_binary_file});
            $self->throw("snp_binary_file must be an absolute path") unless $snp_file->is_absolute;
            my $bin2hapmap = $options->{bin2hapmap_exe};
            my $samtools = $options->{samtools_exe};
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
            	my $tmp_snp_1 = $self->output_file(output_key => 'tmp_files', basename => $bam_basename.'.tmp_1', type => 'txt', temporary => 1);
            	my $tmp_file_1 = $tmp_snp_1->path;
            	my $cmd = qq[$bin2hapmap -l $snp_file > $tmp_file_1 && $samtools mpileup -ugDI -d 1000 -l $tmp_file_1 -f $ref $bam_path > $bcf_path];
                $self->dispatch([$cmd, $req, {output_files => [$bcf_file, $tmp_snp_1]}]);
            }
        return 1;
        };
    }
    
    method outputs_definition {
        return { bcf_files_with_metadata => VRPipe::StepIODefinition->get(type => 'bin',
                                                                  max_files => -1,
                                                                  description => 'bcf file produced for bam file using samtools',
                                                                  metadata => {sample => 'sample name'}
                                                                  ),
                 tmp_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                                  max_files => -1,
                                                                  description => 'tmp files used for samtools and mpileup ',
                                                                  )	
         };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "This step creates a bcf file from a bam file for wgs studies";
    }
    method max_simultaneous {
        return 0;
    }
}

1;
