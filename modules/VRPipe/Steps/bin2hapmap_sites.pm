use VRPipe::Base;

class VRPipe::Steps::bin2hapmap_sites with VRPipe::StepRole {
    method options_definition {
        return { hapmap2bin_sample_genotypes_file => VRPipe::StepOption->get(description => 'absolute path to binary file of sample genotypes, produced by hapmap2bin'),
                 bin2hapmap_exe => VRPipe::StepOption->get(description => 'path to bin2hapmap executable',
                                                           optional => 1,
                                                           default_value => 'bin2hapmap') };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $gtype_snps_bin = Path::Class::File->new($options->{hapmap2bin_sample_genotypes_file});
            $self->throw("hapmap2bin_sample_genotypes_file must be an absolute path") unless $gtype_snps_bin->is_absolute;
            my $bin2hapmap = $options->{bin2hapmap_exe};
            my $req = $self->new_requirements(memory => 3900, time => 1);
            my $hapmap_file = $self->output_file(output_key => 'hapmap_file',
            					 output_dir => $gtype_snps_bin->dir->stringify,
                                                 basename => $gtype_snps_bin->basename.'.hapmap.txt',
                                                 type => 'txt');
            my $hapmap_path = $hapmap_file->path;
            my $cmd = qq[$bin2hapmap -l $gtype_snps_bin > $hapmap_path];
            $self->dispatch([$cmd, $req, {output_files => [$hapmap_file]}]);
        };
    }
    method outputs_definition {
        return { hapmap_file => VRPipe::StepIODefinition->get(type => 'txt',
                                                              description => 'text file with the hapmap sites output from the snp binary file')	};
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Convert a hapmap2bin-generated snp binary file to a textual list of sites with bin2hapmap.";
    }
    method max_simultaneous {
        return 0;
    }
}

1;
