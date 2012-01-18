use VRPipe::Base;
use Data::Dumper;

class VRPipe::Steps::qc_plots_wgs with VRPipe::StepRole {

    method options_definition {
        return { plot_bamcheck_exe => VRPipe::StepOption->get(description => 'path to your plot-bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'plot-bamcheck -p'),
                 reference_fasta_stats => VRPipe::StepOption->get(description => 'absolute path to genome reference stats file used in the plotting process'),
               };
    }
    
    method inputs_definition {
        return { bamcheck_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'bamcheck output files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $stats = Path::Class::File->new($options->{reference_fasta_stats});
            $self->throw("reference_fasta_stats must be an absolute path") unless $stats->is_absolute;
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $plot_bc_exe = $options->{plot_bamcheck_exe};
            foreach my $bc_file (@{$self->inputs->{bamcheck_files}}) {
                my $meta = $bc_file->metadata;
                my $bc_path = $bc_file->path;
                my $graph_dir = $bc_path->dir."/$meta->{sample}";
                my $cmd = qq[$plot_bc_exe $graph_dir $stats $bc_path];
                $self->dispatch([$cmd, $req]);
            }
        };
    }

    method outputs_definition {
        return { };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Uses plot-bamcheck to daw plots from the associated bamcheck file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
