use VRPipe::Base;

class VRPipe::Steps::bam_metadata with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->get(description => 'path to your bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'bamcheck') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', description => 'bam files', max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            my $bamcheck_exe = $options->{bamcheck_exe};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                # our output file is our input file
                my $ifile = $bam_file->path;
                $self->output_file(output_key => 'bam_files_with_metadata', output_dir => $ifile->dir, basename => $ifile->basename, type => 'bam');
                
                my $meta = $bam_file->metadata;
                unless ($meta->{bases} && $meta->{reads} && $meta->{avg_read_length}) {
                    my $check_file = $self->output_file(output_key => 'bamcheck_files', basename => $ifile->basename.'.bamcheck', type => 'txt', temporary => 1);
                    my $ofile = $check_file->path;
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_metadata', 'stats_from_bamcheck', ["$bamcheck_exe $ifile > $ofile", $req, {output_files => [$check_file]}]);
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
                                                                                       reads => 'total number of reads (sequences)',
                                                                                       forward_reads => 'number of forward reads',
                                                                                       reverse_reads => 'number of reverse reads',
                                                                                       avg_read_length => 'the average length of reads',
                                                                                       paired => '0=single ended reads only; 1=paired end reads present',
                                                                                       insert_size => 'average insert size (0 if unpaired)',
                                                                                       library => 'library name',
                                                                                       sample => 'sample name',
                                                                                       center_name => 'center name',
                                                                                       platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                                       study => 'name of the study, put in the DS field of the RG header line',
                                                                                       optional => ['library', 'sample', 'center_name', 'platform', 'study', 'insert_size']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Takes a bam file and associates metadata with the file in the VRPipe database, making the bam file usable in other bam-related Steps";
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
            
            # parse the bamcheck file
            my $parser = VRPipe::Parser->create('bamcheck', {file => $check_file});
            open(my $ifh, $check_path) || die "could not open $check_path\n";
            $new_meta->{reads} = $parser->sequences;
            $new_meta->{paired} = $parser->is_paired;
            $new_meta->{bases} = $parser->total_length;
            $new_meta->{forward_reads} = $parser->first_fragments;
            $new_meta->{reverse_reads} = $parser->last_fragments;
            $new_meta->{avg_read_length} = $parser->average_length;
            $new_meta->{insert_size} = $parser->insert_size_average;
            
            # and get other metadata from bam header
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
                # call the name something we can be most sure is maximally
                # unique
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
