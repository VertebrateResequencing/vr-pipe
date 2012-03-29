use VRPipe::Base;

class VRPipe::Steps::glf_check_genotype with VRPipe::StepRole {
    method options_definition {
        return { hapmap2bin_sample_genotypes_file => VRPipe::StepOption->get(description => 'absolute path to binary file of sample genotypes, produced by hapmap2bin'),
        	 glf_exe => VRPipe::StepOption->get(description => 'path to the glf executable',
                                                    optional => 1,
                                                    default_value => 'glf') };
    }
    method inputs_definition {
        return { bcf_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                            max_files => -1,
                                                            description => 'bcf files for genotyping',
                                                            metadata => {sample => 'name of expected sample', source_bam => 'input bam path'}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $gtype_snps_bin = Path::Class::File->new($options->{hapmap2bin_sample_genotypes_file});
            $self->throw("hapmap2bin_sample_genotypes_file must be an absolute path") unless $gtype_snps_bin->is_absolute;
            my $glf_exe = $options->{glf_exe};
	    
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $bcf (@{$self->inputs->{bcf_files}}) {
		my $bcf_path = $bcf->path;
		my $meta = $bcf->metadata;
		my $sample = $meta->{sample};
		my $source_bam = $meta->{source_bam};
		my $gtypex_file = $self->output_file(output_key => 'gtypex_files_with_metadata',
						     basename => $bcf->basename.'.gtypex',
						     type => 'txt',
						     metadata => {sample => $sample, source_bam => $source_bam});
		my $gtypex_path = $gtypex_file->path;
		my $cmd = qq[$glf_exe checkGenotype -s - $gtype_snps_bin $bcf_path > $gtypex_path];
		$self->dispatch([$cmd, $req, {output_files => [$gtypex_file]}]);	
            }
	};
    } 
    method outputs_definition {
        return { gtypex_files_with_metadata => VRPipe::StepIODefinition->get(type => 'txt',
									     max_files => -1,
									     description => 'file of likelihood scores calculated by glf',
									     metadata => {sample => 'name of expected sample', source_bam => 'input bam path'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Produces gtypex files of genotype likelihood scores for sample genotypes which are calculated by running glf checkGenotype against a snp binary file (of all samples) on bcf files";
    }
    method max_simultaneous {
        return 0;
    }
}

1;
