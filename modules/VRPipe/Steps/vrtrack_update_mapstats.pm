use VRPipe::Base;

class VRPipe::Steps::vrtrack_update_mapstats extends VRPipe::Steps::vrtrack_update {
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                            description => 'bam file with associated bamcheck statistics in the metadata', 
                                                            max_files => -1,
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
									 bases => 'total number of base pairs',
                                                                         reads => 'total number of reads (sequences)',
                                                                         mean_insert_size => 'mean insert size (0 if unpaired)',
                                                                         reads_mapped => 'number of reads mapped',
                                                                         reads_paired => 'number of reads paired',
                                                                         bases_trimmed => 'number of bases trimmed',
                                                                         bases_mapped => 'number of bases mapped',
                                                                         bases_mapped_c => 'number of bases mapped (cigar)',
                                                                         error_rate => 'error rate from bamcheck',
                                                                         rmdup_bases_mapped => 'total number of bases mapped after removing duplicates',
                                                                         rmdup_reads_mapped => 'number of reads mapped after removing duplicates',
                                                                         sd_insert_size => 'standard deviation of the insert size'}),
		 bamcheck_plots => VRPipe::StepIODefinition->get(type => 'bin',
                                                                 description => 'png files produced by plot-bamcheck, with a caption in the metadata',
								 min_files => 11,
                                                                 max_files => -1,
								 metadata => {source_bam => 'the bam file this plot was made from',
									      caption => 'the caption of this plot'}) };
    }
    method body_sub {
        return sub { 
            my $self = shift;
	    my $db = $self->options->{vrtrack_db};
            my $req = $self->new_requirements(memory => 500, time => 1);
	    
	    my %bam_plots;
	    foreach my $plot_file (@{$self->inputs->{bamcheck_plots}}) {
		my $source_bam = $plot_file->metadata->{source_bam};
		$bam_plots{$source_bam}->{dir} = $plot_file->dir;
		push(@{$bam_plots{$source_bam}->{files}}, $plot_file->basename);
	    }
	    
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam_file->path;
		my $lane = $bam_file->metadata->{lane};
		my $these_bam_plots = $bam_plots{$bam_path};
                my $cmd = "use VRPipe::Steps::vrtrack_update_mapstats; VRPipe::Steps::vrtrack_update_mapstats->update_mapstats(db => q[$db], bam => q[$bam_path], lane => q[$lane], plot_dir => q[$these_bam_plots->{dir}], plots => [qw[@{$these_bam_plots->{files}}]]);";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    method description {
        return "Add the bamcheck QC statistics and graphs to the VRTrack database, so that they're accessible with QCGrind etc.";
    }

    method update_mapstats (ClassName|Object $self: Str :$db!, Str|File :$bam!, Str :$lane!, Str|Dir :$plot_dir!, ArrayRef :$plots!) {
	my $bam_file = VRPipe::File->get(path => $bam);  	 
	my $meta = $bam_file->metadata;
	my %plot_files;
	foreach my $plot_basename (@$plots) {
	    my $file = VRPipe::File->get(path => file($plot_dir, $plot_basename));
	    $plot_files{$file->path} = $file->metadata->{caption};
	}
	$bam_file->disconnect;
	
	# do we need to sanity check the .dict file vs the SQ lines in the bam
	# header?
	
	# get the lane and mapstats object from VRTrack
	my $vrtrack = $self->get_vrtrack(db => $db);
	my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("No lane named '$lane' in database '$db'");
	my $mapstats = $vrlane->latest_mapping;
	$vrtrack->transaction_start();
	unless ($mapstats) {
	    $mapstats = $vrlane->add_mapping();
	    #*** do we need to fill in mapper and assembly?
	}
	
	# fill in the mapstats based on our bam metadata
	$mapstats->raw_reads($meta->{reads});
	$mapstats->raw_bases($meta->{bases});
	$mapstats->reads_mapped($meta->{reads_mapped});
	$mapstats->reads_paired($meta->{reads_paired});
	$mapstats->bases_mapped($meta->{bases_mapped_c});
	$mapstats->error_rate($meta->{error_rate});
	$mapstats->rmdup_reads_mapped($meta->{rmdup_reads_mapped});
	$mapstats->rmdup_bases_mapped($meta->{rmdup_bases_mapped}); # ideally this would be rmdup_bases_mapped_c, but we no longer calculate this with bamcheck -d
	$mapstats->clip_bases($meta->{bases} - $meta->{bases_trimmed});
	$mapstats->mean_insert($meta->{mean_insert_size});
	$mapstats->sd_insert($meta->{sd_insert_size});
	
	# add the images
	while (my ($path, $caption) = each %plot_files) {
	    my $img = $mapstats->add_image_by_filename($path);
	    $img->caption($caption);
	    $img->update;
	}
	
	$mapstats->update;
	
	# say that the file is imported
	my $vrfile = $vrlane->get_file_by_name($bam_file->basename);
	if ($vrfile) {
	    $vrfile->is_processed(import => 1);
	    $vrfile->md5($bam_file->md5);
	    $vrfile->update;
	    $vrfile->is_processed(mapped => 1);
	    $vrfile->update;
	}
	
	# also update the lane
	$vrlane->is_processed(import => 1); # we must be imported, but we haven't yet done auto-qc (no genotype), so NOT qc => 1
	$vrlane->raw_bases($meta->{bases});
	$vrlane->raw_reads($meta->{reads});
	$vrlane->is_paired($meta->{paired} ? 1 : 0);
	$vrlane->read_len(int($meta->{avg_read_length}));
	$vrlane->update;
	$vrlane->is_processed(mapped => 1);
	$vrlane->update;
	$vrlane->qc_status('pending');
	$vrtrack->transaction_commit();
    }
}

1;
