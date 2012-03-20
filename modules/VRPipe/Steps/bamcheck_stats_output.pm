use VRPipe::Base;

class VRPipe::Steps::bamcheck_stats_output with VRPipe::StepRole {
    
    method options_definition {
        return { };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', 
                                                            description => 'bam file with associated bamcheck statistics', 
                                                            max_files => -1,
                                                            metadata => {bases => 'total number of base pairs',
                                                                         reads => 'total number of reads (sequences)',
                                                                         mean_insert_size => 'mean insert size (0 if unpaired)',
                                                                         reads_mapped => 'number of reads mapped',
                                                                         reads_paired => 'number of reads paired',
                                                                         bases_trimmed => 'number of bases trimmed',
                                                                         bases_mapped => 'number of bases mapped',
                                                                         bases_mapped_c => 'number of bases mapped (cigar)',
                                                                         error_rate => 'error rate from bamcheck',
                                                                         rmdup_bases => 'number of bases after removing duplicates',
                                                                         rmdup_reads => 'total number of reads after removing duplicates',
                                                                         rmdup_reads_mapped => 'number of reads mapped after removing duplicates',
                                                                         rmdup_bases_mapped_c => 'number of bases mapped after removing duplicates according to the cigar',
                                                                         rmdup_bases_trimmed => 'number of bases trimmed after removing duplicates',
                                                                         sd_insert_size => 'the standard deviation of insert size',
									 optional => ['rmdup_bases_mapped_c', 'rmdup_bases_trimmed']})
               }; 
    }
    
    method body_sub {
        return sub { 
            my $self = shift;
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $stats = $self->output_file(output_key => 'stats_files', basename => $bam_file->basename.'.detailed_stats', type => 'txt');
                my $stats_path = $stats->path;
                my $bam_path = $bam_file->path;
                my $cmd = "use VRPipe::Steps::bamcheck_stats_output; VRPipe::Steps::bamcheck_stats_output->write_stats_file(bam => q[$bam_path], stats => q[$stats_path]);";
                $self->dispatch_vrpipecode($cmd, $req, {output_files => [$stats]});
            }
        };
    }
    
    method outputs_definition {
    	return { stats_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'file to store output from two separate bamcheck statistics, one using the ref fasta and the other removing duplicates', max_files => -1)};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "writes the bamcheck statistics from bam file metadata into a detailed statistics file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }

    method write_stats_file (ClassName|Object $self: Str|File :$bam, Str|File :$stats) {
	my $check_file = VRPipe::File->get(path => $bam);
	my $stats_out = VRPipe::File->get(path => $stats);    	 
	my $meta = $check_file->metadata;
	my $ofh = $stats_out->openw;
	$stats_out->disconnect;
	
	printf $ofh "reads total .. %d\n", $meta->{reads};
	printf $ofh "     mapped .. %d (%.1f%%)\n", $meta->{reads_mapped}, $meta->{reads_mapped}*100./$meta->{reads};
	printf $ofh "     paired .. %d (%.1f%%)\n", $meta->{reads_paired}, $meta->{reads_paired}*100./$meta->{reads};
	printf $ofh "bases total .. %d\n", $meta->{bases};
	printf $ofh "    clip bases     .. %d (%.1f%%)\n", ($meta->{bases}-$meta->{bases_trimmed}),($meta->{bases}-$meta->{bases_trimmed})*100./$meta->{bases};
	printf $ofh "    mapped (read)  .. %d (%.1f%%)\n", $meta->{bases_mapped}, $meta->{bases_mapped}*100./$meta->{bases};
	printf $ofh "    mapped (cigar) .. %d (%.1f%%)\n", $meta->{bases_mapped_c}, $meta->{bases_mapped_c}*100./($meta->{bases}-$meta->{bases_trimmed});
	printf $ofh "error rate  .. %f\n", $meta->{error_rate};
	printf $ofh "rmdup\n";
	printf $ofh "     reads total  .. %d (%.1f%%)\n", $meta->{rmdup_reads}, $meta->{rmdup_reads}*100./$meta->{reads};
	printf $ofh "     reads mapped .. %d (%.1f%%)\n", $meta->{rmdup_reads_mapped}, $meta->{rmdup_reads_mapped}*100./$meta->{rmdup_reads};
	if ($meta->{rmdup_bases_mapped_c} && $meta->{rmdup_bases_trimmed}) {
	   printf $ofh "     bases mapped (cigar) .. %d (%.1f%%)\n",  $meta->{rmdup_bases_mapped_c},  $meta->{rmdup_bases_mapped_c}*100./($meta->{rmdup_bases}-$meta->{rmdup_bases_trimmed});
	}
	printf $ofh "duplication .. %f\n", 1-$meta->{rmdup_reads_mapped}/$meta->{reads_mapped};
	printf $ofh "\n";
	printf $ofh "insert size        \n";
	printf $ofh "    average .. %.1f\n", $meta->{mean_insert_size};
	printf $ofh "    std dev .. %.1f\n", $meta->{sd_insert_size};
	printf $ofh "\n";
        $stats_out->close;
    }
}

1;
