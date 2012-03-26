use VRPipe::Base;
use VRPipe::Parser;

class VRPipe::Steps::bam_snp_sites with VRPipe::StepRole {
	method options_definition {
        return {reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
        		snp_coordinates_file => VRPipe::StepOption->get(description => 'absolute path to file containing the coordinates of the snps locations used in genotyping (samtools will extract these alignments from the bam file provided)'),
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
            my $req = $self->new_requirements(memory => 3900, time => 1);
            my $snp_coordinates_file = Path::Class::File->new($options->{snp_coordinates_file});
            $self->throw("snp_coordinates_file must be an absolute path") unless $snp_coordinates_file->is_absolute;
            my $snp_coordinates_string;
            my $fh = $snp_coordinates_file->openr;
            while (<$fh>) {
                chomp;
                my @a = split /\t/;
                $snp_coordinates_string .= "$a[0]:$a[1]-$a[1] ";
            }
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $bam_basename = $bam->basename;
                $bam_basename =~ s/(.)\.[^.]+$/$1/;
                my $sample;
                my $meta = $bam->metadata;
                if ($meta->{sample}) {
                	$sample = $meta->{sample};
                }
                else {
                	my $parser = VRPipe::Parser->create('bam', {file => $bam_path});
                	my %rg_info = $parser->readgroup_info();
                	my @rgs = keys %rg_info;
                	if (@rgs == 1) {
                		my $info = $rg_info{$rgs[0]};
                		$sample = $info->{SM} if $info->{SM};
                	}
                }
                my $bam_snp_file = $self->output_file(output_key => 'bam_files_with_metadata',
                                              basename => $bam_basename.'.snp.bam',
                                              type => 'bin',
                                              metadata => {sample => $sample});
                my $bam_snp_path = $bam_snp_file->path;
                my $bai_file = $self->output_file(output_key => 'bai_files', basename => $bam_basename.'.snp.bam.bai', type => 'bin');
                my $bai_path = $bai_file->path;
                my $tmp_bam = $self->output_file(basename => $bam_basename.'.snp', type => 'bin', temporary => 1);
                my $tmp_file = $tmp_bam->path;
            	my $cmd = qq[$samtools view -bh $bam_path $snp_coordinates_string > $tmp_file && $samtools sort $tmp_file $tmp_file && $samtools index $bam_snp_path $bai_path];
				$self->dispatch([$cmd, $req, {output_files => [$bam_snp_file, $tmp_bam]}]);
            }
        return 1;
        };
    }
    method outputs_definition {
        return { bam_files_with_metadata => VRPipe::StepIODefinition->get(type => 'bin',
                                                                  max_files => -1,
                                                                  description => 'bam file filtered using snp coordinates file with sample name metadata',
                                                                  metadata => {sample => 'sample name'}
                                                                  ),
                 bai_files => VRPipe::StepIODefinition->get(type => 'bin', 
                                                            max_files => -1, 
                                                            description => 'a .bai file for each bam file')
         };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Bam file is filtered by samtools view using a snp locations file to generate a sorted and indexed bam file of sample bam snp sites";
    }
    method max_simultaneous {
        return 0;
    }
}
1;
