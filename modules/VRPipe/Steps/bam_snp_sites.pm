use VRPipe::Base;

class VRPipe::Steps::bam_snp_sites with VRPipe::StepRole {

	method options_definition {
        return {reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
        		snp_sites_file => VRPipe::StepOption->get(description => 'absolute path to snp sites file used for genotyping'),
			    samtools_view_options => VRPipe::StepOption->get(description => 'samtools view options', 
					optional => 1, 
					default_value => '-bh'),
        		samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools')
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
            my $samtools = $options->{samtools_exe};
            my $samtools_view_opts = $options->{samtools_view_options};
            my $req = $self->new_requirements(memory => 3900, time => 1);
            my $snp_file = Path::Class::File->new($options->{snp_sites_file});
            $self->throw("snp_sites_file must be an absolute path") unless $snp_file->is_absolute;
            my $snp_sites_string;
            $self->throw("Error opening file ". $snp_file) unless open (SNPFILE, "<", $snp_file);
            while (<SNPFILE>) {
                chomp;
                my @a = split /\t/;
                $snp_sites_string .= "$a[0]:$a[1]-$a[1] ";
            }
            close SNPFILE;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $bam_dir = $bam->dir;
                my $bam_basename = $bam->basename;
                $bam_basename =~ s/(.)\.[^.]+$/$1/;
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
                my $bam_snp_file = $self->output_file(output_key => 'bam_files_with_metadata',
                                              basename => $bam_basename.'.snp.bam',
                                              type => 'bin',
                                              metadata => {sample => $sample});
                my $bam_snp_path = $bam_snp_file->path;
                my $tmp_bam = $self->output_file(output_key => 'tmp_files', basename => $bam_basename.'.snp', type => 'bin', temporary => 1);
            	my $tmp_file = $tmp_bam->path;
            	my $cmd = qq[$samtools view $samtools_view_opts $bam_path $snp_sites_string > $tmp_file && $samtools sort $tmp_file $tmp_file && $samtools index $bam_snp_path];
				$self->dispatch([$cmd, $req, {output_files => [$bam_snp_file, $tmp_bam]}]);
            }
        return 1;
        };
    }
    method outputs_definition {
        return { bam_files_with_metadata => VRPipe::StepIODefinition->get(type => 'bin',
                                                                  max_files => -1,
                                                                  description => 'bam file filtered using snp sites file with sample name metadata',
                                                                  metadata => {sample => 'sample name'}
                                                                  ),
                 tmp_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                                  max_files => -1,
                                                                  description => 'tmp files used for samtools and mpileup ',
                                                                  )		
         };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Bam file is filtered by samtools view using a snp sites file to generate a sorted and indexed bam file of sample bam snp sites";
    }
    method max_simultaneous {
        return 0;
    }
}
1;
