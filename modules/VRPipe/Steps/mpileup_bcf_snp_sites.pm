use VRPipe::Base;

class VRPipe::Steps::mpileup_bcf_snp_sites extends VRPipe::Steps::mpileup_bcf {
    method inputs_definition {
        return { 
        		 bam_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                            max_files => -1,
                                                            description => 'bam files with sample name metadata',
                                                            metadata => {sample => 'sample name'}
                                                            )
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
            my $snp_file = Path::Class::File->new($options->{snp_sites_file});
            $self->throw("snp_sites_file must be an absolute path") unless $snp_file->is_absolute;
            my $mpileup_opts = $options->{samtools_mpileup_options};
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $bam_dir = $bam->dir;
                my $bam_basename = $bam->basename;
                my $meta = $bam->metadata;
                my $sample = $meta->{sample};
                my $bcf_file = $self->output_file(output_key => 'bcf_files_with_metadata',
                                              basename => $bam_basename.'.bcf',
                                              type => 'bin',
                                              metadata => {sample => $sample});
                my $bcf_path = $bcf_file->path;
            	my $cmd = qq[$samtools mpileup -$mpileup_opts $snp_file -f $ref $bam_path > $bcf_path];
				$self->dispatch([$cmd, $req, {output_files => [$bcf_file]}]);#, $tmp_index]}]);
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
         };
    }
    method description {
        return "Run samtools mpileup for each bam using the snps sites file provided to generate one bcf file with associated sample name metadata per bam";
    }
}
1;
