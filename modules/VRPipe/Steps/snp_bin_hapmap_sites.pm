use VRPipe::Base;

class VRPipe::Steps::snp_bin_hapmap_sites with VRPipe::StepRole {

	method options_definition {
        return { snp_binary_file => VRPipe::StepOption->get(description => 'absolute path to snp binary                                                    file used for genotyping'),
        		 bin2hapmap_exe => VRPipe::StepOption->get(description => 'path to bin2hapmap                                                executable',
                                                           optional => 1,
                                                           default_value => 'bin2hapmap'),
                 bin2hapmap_options => VRPipe::StepOption->get(description => 'options provided for                                               bin2hapmap',
                                                           optional => 1,
                                                           default_value => '-l')
        };
    }

    method inputs_definition {
        return { 
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $snp_file = Path::Class::File->new($options->{snp_binary_file});
            $self->throw("snp_binary_file must be an absolute path") unless $snp_file->is_absolute;
            my $snp_basename = $snp_file->basename;
            my $bin2hapmap = $options->{bin2hapmap_exe};
            my $bin2hapmap_options = $options->{bin2hapmap_options};
            my $req = $self->new_requirements(memory => 3900, time => 1);
            my $hapmap_file = $self->output_file(output_key => 'hapmap_file',
            									 output_dir => $snp_file->dir->stringify,
                                                 basename => $snp_basename.'.hapmap.txt',
                                                 type => 'txt');
            my $hapmap_path = $hapmap_file->path;
            my $cmd = qq[$bin2hapmap $bin2hapmap_options $snp_file > $hapmap_path];
            $self->dispatch([$cmd, $req, {output_files => [$hapmap_file]}]);
            return 1;
        };
    }
    
    method outputs_definition {
        return { hapmap_file => VRPipe::StepIODefinition->get(type => 'txt',
                                                              description => 'text file with the hapmap sites output from the snp binary file')	
         };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "This step uses bin2hapmap to create a hapmap sites file from a snp binary file";
    }
    method max_simultaneous {
        return 0;
    }
}

1;
