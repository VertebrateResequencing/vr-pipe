use VRPipe::Base;

class VRPipe::Steps::genotype_checking with VRPipe::StepRole {

	method options_definition {
        return { snp_binary_file => VRPipe::StepOption->get(description => 'absolute path to binary file used for genotyping'),
        		 glf_exe => VRPipe::StepOption->get(description => 'path to your glf executable',
                                                    optional => 1,
                                                    default_value => 'glf')
               };
    }
    
    method inputs_definition {
        return { 
        		 bcf_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                            max_files => -1,
                                                            description => 'bcf files for genotyping',
                                                            metadata => {sample => 'sample name'}
                                                            )
        		 };
	}
	
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $snp_bin = Path::Class::File->new($options->{snp_binary_file});
            $self->throw("snp_binary_file must be an absolute path") unless $snp_bin->is_absolute;
            my $glf_exe = $options->{glf_exe};
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $bcf (@{$self->inputs->{bcf_files}}) {
                my $bcf_path = $bcf->path;
                my $bcf_basename = $bcf->basename;
                my $meta = $bcf->metadata;
                my $sample;
                if ($meta->{sample}) {
                	$sample = $meta->{sample};
                }
                $self->throw("unable to obtain sample name metadata from the bcf file") unless $sample;
                my $gtypex_file = $self->output_file(output_key => 'gtypex_files_with_metadata',
                                              basename => $bcf_basename.'.gtypex',
                                              type => 'txt',
                                              metadata => {sample => $sample});
                my $gtypex_path = $gtypex_file->path;
            	my $cmd = qq[$glf_exe checkGenotype -s - $snp_bin $bcf_path > $gtypex_path];
                $self->dispatch([$cmd, $req, {output_files => [$gtypex_file]}]);	
            }
       };
   } 

    method outputs_definition {
        return { gtypex_files_with_metadata => VRPipe::StepIODefinition->get(type => 'txt',
                                                                  max_files => -1,
                                                                  description => 'file of likelihood scores calculated by glf',
                                                                  metadata => {sample => 'sample name'}
                                                                  )	
         };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "This step performs checkGenotype using glf on bcf files";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
