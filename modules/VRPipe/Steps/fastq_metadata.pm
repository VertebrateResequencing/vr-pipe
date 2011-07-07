use VRPipe::Base;

class VRPipe::Steps::fastq_metadata with VRPipe::StepRole {
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq', description => 'fastq files', max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            foreach my $fq_file (@{$self->inputs->{fastq_files}}) {
                # our output file is our input file
                my $ifile = $fq_file->path;
                $self->output_file(output_key => 'fastq_files', output_dir => $ifile->dir, basename => $ifile->basename);
                
                # we expect the datasource to fill in most if not all metadata,
                # but some things may need calculating
                my $meta = $fq_file->metadata;
                unless ($meta->{bases} && $meta->{reads}) {
                    # run fastqcheck to generate stats
                    # record the fact that we (probably) made the fastqcheck
                    # file, even though it isn't in our outputs_definition *** or should we be able to have optional outputs?
                    my $fqc_file = $self->output_file(output_key => 'fastqcheck_files', output_dir => $ifile->dir, basename => $ifile->basename.'.fastqcheck', type => 'txt', temporary => 1);
                    my $ofile = $fqc_file->path;
                    unless ($fqc_file->s) {
                        my $req = $self->new_requirements(memory => 50, time => 1);
                        if ($ifile =~ /\.gz/) {
                            $self->dispatch([qq{gunzip -c $ifile | fastqcheck > $ofile}, $req]);
                        }
                        else {
                            $self->dispatch([qq{fastqcheck $ifile > $ofile}, $req]);
                        }
                    }
                }
                unless ($fq_file->md5) {
                    $self->dispatch_md5sum($fq_file, $meta->{expected_md5});
                }
            }
        };
    }
    method outputs_definition {
        return { fastq_files_with_metadata => VRPipe::StepIODefinition->get(type => 'fq',
                                                                            description => 'a fastq file with associated metadata',
                                                                            max_files => -1,
                                                                            metadata => {bases => 'total number of base pairs',
                                                                                         reads => 'total number of reads (sequences)',
                                                                                         #*** etc.
                                                                                         }) };
    }
    method post_process_sub {
        return sub {
            my $self = shift;
            
            my $all_ok = 1;
            foreach my $ofile (@{$self->outputs->{fastq_files}}) {
                my $meta = $ofile->metadata;
                unless ($meta->{bases} && $meta->{reads}) {
                    # body_sub should have created a .fastqcheck file
                    my $fqc_file = VRPipe::File->get(path => $ofile->path.'.fastqcheck');
                    if ($fqc_file->s) {
                        # parse the file
                        my $parser = VRPipe::Parser->create('fqc', {file => $fqc_file});
                        $meta->{bases} = $parser->total_length;
                        $meta->{reads} = $parser->num_sequences;
                        #*** etc.
                        
                        $ofile->add_metadata($meta);
                    }
                    else {
                       $self->warn("some metadata was missing from the dataelement of ".$ofile->path.", and so is the .fastqcheck that could have resolved that");
                       $all_ok = 0;
                    }
                }
                
                unless ($ofile->md5) {
                    $self->warn("missing md5 for ".$ofile->path);
                    $all_ok = 0;
                }
            }
            
            return $all_ok;
        };
    }
    method description {
        return "Takes a fastq file and associates metadata with the file in the VRPipe database, making the fastq file usable in other fastq-related Steps";
    }
}

1;
