use VRPipe::Base;
use Data::Dumper;
class VRPipe::Steps::qc_stats_bamcheck_wgs with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->get(description => 'path to your bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'bamcheck -q 20 -r'),
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'), };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', description => 'bam files', max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift; 
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $bamcheck_exe = $options->{bamcheck_exe};
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $ifile = $bam_file->path;
                $self->output_file(output_key => 'bam_files_with_metadata', output_dir => $ifile->dir, basename => $ifile->basename, type => 'bam');
                my $meta = $bam_file->metadata; 
                unless ($meta->{bases} && $meta->{reads} && $meta->{avg_read_length}) {
                    my $check_file = $self->output_file(output_key => 'bamcheck_files', basename => $ifile->basename.'.bamcheck', type => 'txt', temporary => 1);
                    my $ofile = $check_file->path;
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::qc_stats_bamcheck_wgs', 'stats_from_bamcheck', ["$bamcheck_exe $ref $ifile > $ofile", $req, {output_files => [$check_file]}]);
                }
            }
        };
    }
    method outputs_definition {
        return { bam_files_with_metadata => VRPipe::StepIODefinition->get(type => 'bam',
                                                                          description => 'a bam file with associated metadata',
                                                                          max_files => -1,
                                                                          metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                          bases => 'total number of base pairs',
                                                                          bases_mapped => 'number of base pairs mapped',
                                                                          bases_mapped_c => 'number of base pairs mapped (cigar)',
                                                                          bases_trimmed => 'number of trimmed base pairs',
                                                                          reads => 'total number of reads (sequences)',
                                                                          reads_mapped => 'total number of reads (sequences) mapped',
                                                                          reads_paired => 'total number of reads (sequences) paired',
                                                                          insert_size => 'the average insert size (0 if unpaired)',
                                                                          sd_insert_size => 'sd of average insert size',
                                                                          error_rate => 'error rate of sequencing',
                                                                          library => 'library name',
                                                                          sample => 'sample name',
                                                                          center_name => 'center name',
                                                                          platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                          study => 'name of the study, put in the DS field of the RG header line',
                                                                          optional => ['library', 'sample', 'center_name', 'platform', 'study']}) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Takes a bam file and associates metadata with the file for production of stats and graphs by later steps";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method stats_from_bamcheck (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $check_path) = $cmd_line =~ / (\S+) > (\S+)$/;
        $bam_path || $self->throw("bad cmd line [$cmd_line]");
        my $bam_file = VRPipe::File->get(path => $bam_path);
        my $check_file = VRPipe::File->get(path => $check_path);
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $check_file->update_stats_from_disc(retries => 3);
        if ($check_file->s) {
            my $new_meta = {};
            # parse the bamcheck file to get stats data
            my $parser = VRPipe::Parser->create('bamcheck', {file => $check_file});
            open(my $ifh, $check_path) || die "could not open $check_path\n";
    		$new_meta->{reads} = $parser->sequences;
			$new_meta->{reads_mapped} = $parser->reads_mapped;
			$new_meta->{reads_paired} = $parser->reads_paired;
			$new_meta->{bases} = $parser->total_length;
			$new_meta->{bases_mapped} = $parser->bases_mapped;
			$new_meta->{bases_mapped_c} = $parser->bases_mapped_cigar;
			$new_meta->{bases_trimmed} = $parser->bases_trimmed;
			$new_meta->{error_rate} = $parser->error_rate;
			$new_meta->{insert_size} = $parser->insert_size_average;
			$new_meta->{sd_insert_size} = $parser->insert_size_standard_deviation;
            # get metadata from bam header
            $parser = VRPipe::Parser->create('bam', {file => $bam_file});
            my %rg_info = $parser->readgroup_info();
            my @rgs = keys %rg_info;
            if (@rgs == 1) {
                my $info = $rg_info{$rgs[0]};
                $new_meta->{lane} = $info->{PU} ? $info->{PU} : $rgs[0];
                $new_meta->{library} = $info->{LB} if $info->{LB};
                $new_meta->{sample} = $info->{SM} if $info->{SM};
                $new_meta->{center_name} = $info->{CN} if $info->{CN};
                $new_meta->{platform} = $info->{PL} if $info->{PL};
                $new_meta->{study} = $info->{DS} if $info->{DS};
            }
            if (@rgs != 1 || "$new_meta->{lane}" eq "1") {
                my $rg = $bam_path;
                $rg =~ s/\.bam$//;
                $rg =~ s/\W/_/g;
                $new_meta->{lane} = $rg;
            }
            $bam_file->add_metadata($new_meta);
        }
        else {
           $self->throw("$check_path failed to be made");
        }
    }
}

1;
