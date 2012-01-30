use VRPipe::Base;

class VRPipe::Steps::qc_plots_exome with VRPipe::StepRole {
    use VRPipe::Parser;
    use VRPipe::Parser::fastq;
    use VRPipe::Utils::Math;
    use VRPipe::Utils::Graphs;
	use List::Util qw(min max sum);

	method options_definition {
        return { plot_type => VRPipe::StepOption->get(description => 'type of plot required', optional => 1, default_value => 'png') };
    }
    
    method inputs_definition {
        return { stats_dump_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'dump file to store result of qc_stats_exome', max_files => -1, metadata => {sample => 'sample name'} 
               ) }; 
    }
    
    method body_sub {
        return sub { 
            my $self = shift;
            my $options = $self->options;
            my $plot_type = $options->{plot_type};
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $dump (@{$self->inputs->{stats_dump_files}}) {
                my $meta = $dump->metadata;
                my $out_prefix = $meta->{sample};
                my $dump_path = $dump->path;
                my $cmd = "use VRPipe::Steps::qc_plots_exome; VRPipe::Steps::qc_plots_exome->bam_exome_qc_make_plots(dump => q[$dump_path], outfiles_prefix => q[$out_prefix], plot_type => q[$plot_type]);";
                $self->dispatch_vrpipecode($cmd, $req);
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
        return "Produces plots for sample bam file from a statistics file using VRPipe::Utils::Graphs";
    }
    method max_simultaneous {
        return 0;
    }
    
=head2 bam_exome_qc_make_plots

 Title   : bam_exome_qc_make_plots
 Usage   : $obj->bam_exome_qc_make_plots(\%stats, 'outfiles_prefix', 'plot_type');
 Function: Makes exome QC plots using data gathered by bam_exome_qc_stats
 Returns : nothing (makes plots)
 Args    : reference to hash (returned by bam_exome_qc_stats)
           prefix of output plot files to be made
           plot type (must be recognised by R, e.g. 'pdf', 'png', ...)

=cut
    method bam_exome_qc_make_plots (
     ClassName|Object $self: Str|File :$dump, Str :$outfiles_prefix, Str :$plot_type) {
		my $stats = do $dump or $self->throw("Could not load stats from file $dump");
    	my %coverage;
    	my %plot_data;
        my $math_util = VRPipe::Utils::Math->new();
        my $graph_util = VRPipe::Utils::Graphs->new();
    	foreach my $bait_or_target qw(bait target) {
        	foreach my $hash (@{$stats->{'gc_vs_' . $bait_or_target . '_cvg'}}) {
            	while (my ($cov, $freq) = each(%{$hash})) {
                	$coverage{$bait_or_target}{$cov} += $freq;
            	}
        	}
            my %cov_stats = $math_util->histogram_stats($coverage{$bait_or_target});
            my $xmin = max (0, $cov_stats{mean} - 4 * $cov_stats{standard_deviation});
        	my $xmax = $cov_stats{mean} + 4 * $cov_stats{standard_deviation};
        	$plot_data{$bait_or_target}{x} = [];
        	$plot_data{$bait_or_target}{y} = [];

        	foreach my $i ($xmin .. $xmax) {
            	if (exists $coverage{$bait_or_target}{$i}) {
                	push @{$plot_data{$bait_or_target}{x}}, $i;
                	push @{$plot_data{$bait_or_target}{y}}, $coverage{$bait_or_target}{$i};
            	}
        	}
    	}
    	
    	my @data = ({xvals => $plot_data{bait}{x}, yvals => $plot_data{bait}{y}, legend => 'baits'},
                {xvals => $plot_data{target}{x}, yvals => $plot_data{target}{y}, legend => 'targets'},);
        
        my $msg = $graph_util->plot_stats({outfile=>"$outfiles_prefix.mean_coverage.$plot_type",
                        title => "Mean coverage distribution of baits and targets",
                        desc_xvals => 'Mean coverage',
                        desc_yvals => 'Frequency',
                        data => \@data});
                        
    	# Coverage plots
    	%plot_data = ();

        foreach my $bait_or_target qw(bait target) {
        	my %cov_stats = $math_util->histogram_stats($stats->{$bait_or_target . '_cvg_hist'});
        	my $xmin = max (0, $cov_stats{mean} - 4 * $cov_stats{standard_deviation});
            my $xmax = $cov_stats{mean} + 4 * $cov_stats{standard_deviation};

        	$plot_data{$bait_or_target}{x_per_base} = [];
        	$plot_data{$bait_or_target}{y_per_base} = [];
        	$plot_data{$bait_or_target}{x_norm} = [];
        	$plot_data{$bait_or_target}{y_norm} = [];
        	$plot_data{$bait_or_target}{x_cumulative} = [];
        	$plot_data{$bait_or_target}{y_cumulative} = [];

        	# get the per-base coverage data
        	foreach my $i ($xmin .. $xmax) {
            	if (exists $stats->{$bait_or_target . '_cvg_hist'}{$i}) {
                	push @{$plot_data{$bait_or_target}{x_per_base}}, $i;
                	push @{$plot_data{$bait_or_target}{y_per_base}}, $stats->{$bait_or_target . '_cvg_hist'}{$i};
            	}
        	}

        	foreach my $cov (sort {$a <=> $b} keys %{$stats->{$bait_or_target . '_cumulative_coverage_pct'}}) {
            	push @{$plot_data{$bait_or_target}{x_cumulative}}, $cov;
            	push @{$plot_data{$bait_or_target}{y_cumulative}}, $stats->{$bait_or_target . '_cumulative_coverage_pct'}{$cov};
            	last if $cov >= 2 * $stats->{mean_target_coverage};
        	}

        	# work out the normalised coverage plot data
        	foreach my $cov (sort {$a <=> $b} keys %{$stats->{$bait_or_target . '_cumulative_coverage'}}) {
            	my $bases_fraction = $stats->{$bait_or_target . '_cumulative_coverage'}{$cov} / $stats->{$bait_or_target . '_bases'};
            	my $normalised_cov = $cov / $stats->{'mean_' . $bait_or_target . '_coverage'};
            	push @{$plot_data{$bait_or_target}{x_norm}}, $normalised_cov;
            	push @{$plot_data{$bait_or_target}{y_norm}}, $bases_fraction;
            	last if $normalised_cov > 1;
        	}
    	}
    	@data = ({xvals => $plot_data{bait}{x_per_base}, yvals => $plot_data{bait}{y_per_base}, legend => 'baits'},
                {xvals => $plot_data{target}{x_per_base}, yvals => $plot_data{target}{y_per_base}, legend => 'targets'},);

    	$graph_util->plot_stats({outfile=>"$outfiles_prefix.coverage_per_base.$plot_type",
                        	title => "Bait and target coverage per base",
                        	desc_xvals => 'Coverage',
                        	desc_yvals => 'Frequency',
                        	data => \@data});

    	@data = ({xvals => $plot_data{bait}{x_norm}, yvals => $plot_data{bait}{y_norm}, legend => 'baits'},
                {xvals => $plot_data{target}{x_norm}, yvals => $plot_data{target}{y_norm}, legend => 'targets'},);

        $graph_util->plot_stats({outfile=>"$outfiles_prefix.normalised_coverage.$plot_type",
                           title => "Bait and target normalised coverage",
                           desc_xvals => 'Normalised coverage',
                           desc_yvals => 'Fraction of bases',
                           r_plot => 'ylim=c(0,1)', 
                           data => \@data});

    	@data = ({xvals => $plot_data{bait}{x_cumulative}, yvals => $plot_data{bait}{y_cumulative}, legend => 'baits'},
                {xvals => $plot_data{target}{x_cumulative}, yvals => $plot_data{target}{y_cumulative}, legend => 'targets'},);

        $graph_util->plot_stats({outfile=>"$outfiles_prefix.cumulative_coverage.$plot_type",
                           title => "Bait and target cumulative coverage",
                           desc_xvals => 'Coverage',
                           desc_yvals => 'Fraction of bases',
                           r_plot => 'ylim=c(0,1)', 
                           data => \@data});

    	# Plot: GC of baits/targets vs mean, quartiles coverage
    	foreach (qw(bait target)) {
        	$graph_util->plot_histograms_distributions({outfile=>"$outfiles_prefix.$_" . "_gc_vs_cvg.$plot_type",
                                                   title => "Coverage vs GC plot of $_" . 's',
                                               	   desc_xvals => 'GC (%)',
                                                   desc_yvals => 'Mapped depth',
                                                   xdata => [(0..100)],
                                                   ydata => $stats->{"gc_vs_$_" . '_cvg'}});
        	$graph_util->plot_histograms_distributions({outfile=>"$outfiles_prefix.$_" . "_gc_vs_cvg.scaled.$plot_type",
                                                   title => "Coverage vs GC plot of $_" . 's',
                                                   x_scale => 1, 
                                                   x_scale_values => [(30,40,50)],
                                                   desc_xvals => "Percentile of $_ sequence ordered by GC content",
                                                   desc_xvals_top => '%GC',
                                                   desc_yvals => 'Mapped depth',
                                                   xdata => [(0..100)],
                                                   ydata => $stats->{"gc_vs_$_" . '_cvg'}});
        }

    	# Plot: Quality scores by cycle
    	foreach (qw(1 2 up)) {
        	my %key2string = (1, 'first of pair', 2, 'second of pair', 'up', 'unpaired');
        	$graph_util->plot_histograms_distributions({outfile=>"$outfiles_prefix.quality_scores_$_.$plot_type",
                                              	   xdata => [(0 .. $stats->{readlen} - 1)],
                                                   ydata => $stats->{'qual_scores_' . $_},
                                                   title => "Quality scores distribution, $key2string{$_}",
                                                   desc_xvals => 'Cycle',
                                                   desc_yvals => 'Quality score',
                                                   y_min => 0,
                                                   y_max => 41});
    	}

    	# Plot: reads/targets/baits GC content
    	my @read_xvals = (0..$stats->{readlen});
    	foreach(@read_xvals){$_ = 100 * $_ / $stats->{readlen}}

    	my @zero_to_100 = (0..100);

    	foreach my $type ('unmapped', 'mapped') {
        	my @yvals_1 = @{$stats->{"gc_hist_$type" . '_1'}};
        	my @yvals_2 = @{$stats->{"gc_hist_$type" . '_2'}};
        	my @data = ({xvals => \@read_xvals, yvals => \@yvals_1, legend => '_1'},
            	        {xvals => \@read_xvals, yvals => \@yvals_2, legend => '_2'},
                	    {xvals => \@zero_to_100, yvals => $stats->{bait_gc_hist}, legend => 'bait'},
                    	{xvals => \@zero_to_100, yvals => $stats->{target_gc_hist}, legend => 'target'},);
        	$graph_util->plot_stats({outfile=>"$outfiles_prefix.gc_$type.$plot_type",
                            	normalize => 1,
                            	title => "GC Plot $type reads",
                            	desc_xvals => '%GC',
                            	desc_yvals => 'Frequency',
                            	data => \@data});
    	}

    	# Plot: insert size
    	my %insert_stats = $math_util->new()->histogram_stats($stats->{insert_size_hist});

    	# plot 4 standard devations away from the mean
    	my $xmin = max (0, $insert_stats{mean} - 4 * $insert_stats{standard_deviation});
    	my $xmax = $insert_stats{mean} + 4 * $insert_stats{standard_deviation};
    	my @xvals = ();
    	my @yvals = ();

    	foreach my $i ($xmin .. $xmax) {
        	if (exists $stats->{insert_size_hist}{$i}) {
	            push @xvals, $i;
    	        push @yvals, $stats->{insert_size_hist}{$i};
        	}
    	}

    	$graph_util->plot_stats({outfile=>"$outfiles_prefix.insert_size.$plot_type",
                        	title => "Insert Size Distribution",
                        	desc_xvals => 'Insert Size',
                        	desc_yvals => 'Number of read pairs',
                        	data => [{xvals => \@xvals, yvals =>\@yvals}]});
	}
}

1;
