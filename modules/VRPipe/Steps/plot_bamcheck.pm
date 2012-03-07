use VRPipe::Base;

class VRPipe::Steps::plot_bamcheck with VRPipe::StepRole {

    method options_definition {
        return { plot_bamcheck_exe => VRPipe::StepOption->get(description => 'path to plot-bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'plot-bamcheck'),
                 plot_bamcheck_options => VRPipe::StepOption->get(description => 'options for plot-bamcheck, excluding -p and -r',
                                                         optional => 1),
                 reference_fasta_stats => VRPipe::StepOption->get(description => 'absolute path to genome reference stats file used in the plotting process', optional => 1),
               };
    }
    
    method inputs_definition {
        return { bamcheck_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'bamcheck output files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $plot_bc_exe = $options->{plot_bamcheck_exe};
            my $plot_bc_opts = $options->{plot_bamcheck_options};
	    if ($plot_bc_opts =~ /-p|-r/) {
		$self->throw("plot_bamcheck_options must not contain -p or -r");
	    }
            if ($options->{reference_fasta_stats}) {
            	my $ref_stats = Path::Class::File->new($options->{reference_fasta_stats});
            	$self->throw("reference_fasta_stats must be an absolute path") unless $ref_stats->is_absolute;
		$plot_bc_opts .= '-r '.$ref_stats;
            }
	    
            foreach my $bc_file (@{$self->inputs->{bamcheck_files}}) {
                my $meta = $bc_file->metadata;
		my $prefix = $meta->{lane} || $bc_file->basename;
		my $source_bam = $meta->{source_bam};
		
		$self->output_file(temporary => 1, basename => $prefix.'-quals.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-quals.png', type => 'bin', metadata => {caption => 'Qualities', source_bam => $source_bam});
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-quals2.png', type => 'bin', metadata => {caption => 'Qualities', source_bam => $source_bam});
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-quals3.png', type => 'bin', metadata => {caption => 'Qualities', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-quals-hm.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-quals-hm.png', type => 'bin', metadata => {caption => 'Qualities', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-insert-size.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-insert-size.png', type => 'bin', metadata => {caption => 'Insert Size', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-gc-content.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-gc-content.png', type => 'bin', metadata => {caption => 'GC Content', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-gc-depth.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-gc-depth.png', type => 'bin', metadata => {caption => 'GC Depth', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-acgt-cycles.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-acgt-cycles.png', type => 'bin', metadata => {caption => 'ACGT Cycles', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-coverage.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-coverage.png', type => 'bin', metadata => {caption => 'Coverage', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-indel-dist.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-indel-dist.png', type => 'bin', metadata => {caption => 'Indel distribution', source_bam => $source_bam});
		$self->output_file(temporary => 1, basename => $prefix.'-indel-cycles.gp', type => 'txt');
		$self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-indel-cycles.png', type => 'bin', metadata => {caption => 'Indels per cycle', source_bam => $source_bam});
		
		#*** really we need to know if the bamcheck file was produced using 'bamcheck -r'
		if ($options->{reference_fasta_stats}) {
		    $self->output_file(temporary => 1, basename => $prefix.'-mism-per-cycle.gp', type => 'txt');
		    $self->output_file(output_key => 'bamcheck_plots', basename => $prefix.'-mism-per-cycle.png', type => 'bin', metadata => {caption => 'Mismatches per cycle', source_bam => $source_bam});
		}
		
		my $bc_path = $bc_file->path;
                my $cmd = qq[$plot_bc_exe $plot_bc_opts -p $prefix $bc_path];
                $self->dispatch([$cmd, $req]);
            }
        };
    }

    method outputs_definition {
        return { bamcheck_plots => VRPipe::StepIODefinition->get(type => 'bin',
                                                                 description => 'png files produced by plot-bamcheck, with a caption in the metadata',
								 min_files => 11,
                                                                 max_files => -1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Uses plot-bamcheck to draw plots from the associated bamcheck file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
