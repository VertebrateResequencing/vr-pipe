use VRPipe::Base;

class VRPipe::Steps::fastq_import with VRPipe::StepRole {
    use File::Fetch;
    use Net::FTP::Robust;
    
    method options_definition {
        return { };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                              description => 'fastq files',
                                                              max_files => -1,
                                                              metadata => {remote_path => 'the complete remote location of the file',
                                                                           expected_md5 => 'the md5 checksum the file is supposed to have',
                                                                           optional => ['remote_path', 'expected_md5']},
                                                              check_existence => 0) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $fq_file (@{$self->inputs->{fastq_files}}) {
                # our output file is our input file
                my $ifile = $fq_file->path;
                $self->output_file(output_key => 'local_fastq_files', output_dir => $ifile->dir, basename => $ifile->basename, type => 'fq');
                
                my $meta = $fq_file->metadata;
                my $remote_path = $meta->{remote_path};
                if ($remote_path && $remote_path ne $ifile) {
                    $self->dispatch_vrpipecode(qq[use VRPipe::Steps::fastq_import; VRPipe::Steps::fastq_import->download_fastq(source => q[$remote_path], dest => q[$ifile]);],
                                               $req);
                }
            }
        };
    }
    method outputs_definition {
        return { local_fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                                    description => 'a fastq file on a local disc',
                                                                    max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "If fastq files in the datasource are on an external ftp site, downloads them to local disc";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method download_fastq (ClassName|Object $self: Str :$source, Str|File :$dest) {
        my $fq_file = VRPipe::File->get(path => $dest);
        my $meta = $fq_file->meta;
        my $expected_md5 = $meta->{expected_md5};
        my $local_dir = $fq_file->dir;
        my $out_file = $fq_file->path;
        
        $fq_file->disconnect;
        
        my $ff = File::Fetch->new(uri => $source);
        my $scheme = $ff->scheme;
        my $host = $ff->host;
        my $path = $ff->path;
        my $basename = $ff->file;
        my $full_path = $path.$basename;
        
        if ($scheme eq 'ftp') {
            # use Net::FTP::Robust, since it's potentially better
            my $ftp = Net::FTP::Robust->new(Host => $host);
            
            $ftp->get($full_path, $local_dir);
            
            $fq_file->update_stats_from_disc;
            unless ($fq_file->s) {
                $self->throw("After 10 automated attempts, failed to download $source at all");
            }
            
            if ($expected_md5) {
                my $ok = $self->verify_md5($out_file, $expected_md5);
                
                unless ($ok) {
                    my $tries = 0;
                    while (! $ok) {
                        $fq_file->unlink;
                        
                        $fq_file->disconnect;
                        $ftp->get($full_path, $local_dir);
                        $fq_file->update_stats_from_disc;
                        $ok = $self->verify_md5($out_file, $expected_md5);
                        
                        $tries++;
                        last if $tries == 3;
                    }
                }
                
                unless ($ok) {
                    $fq_file->unlink;
                    $self->throw("Tried downloading $source 3 times, but the md5 never matched '$expected_md5'");
                }
            }
        }
        else {
            $self->throw("'$scheme' downloads are not yet supported");
        }
    }
}

1;
