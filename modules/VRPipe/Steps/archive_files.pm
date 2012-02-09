use VRPipe::Base;

class VRPipe::Steps::archive_files with VRPipe::StepRole {
    use Digest::MD5;
    
    use vars qw(%pool_dirs @robin $robin_index); # we can't use 'our' because of issues when we eval body_sub
    
    method options_definition {
        return { disc_pool_file => VRPipe::StepOption->get(description => 'path to a file with an absolute path of a storage root directory on each line')};
    }
    method inputs_definition {
        return { files => VRPipe::StepIODefinition->get(type => 'any', description => 'files', max_files => -1) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            my $disc_pool_file = Path::Class::File->new($options->{disc_pool_file});
            $self->throw("disc_pool_file must be an absolute path") unless $disc_pool_file->is_absolute;
            
            # parse the disc_pool_file; we don't cache these results in a class
            # variable because disc_pool_file could be updated at any time.
            # we do keep track of all dirs on the class so we maintain a global
            # pool of output locations we can round-robin
            my $fh = $disc_pool_file->openr;
            my %dirs;
            while (<$fh>) {
                chomp;
                my $dir = $_;
                next unless dir($dir)->is_absolute;
                $dirs{$dir} = 1;
                $pool_dirs{$dir} = 1;
            }
            
            # update robin
            my %current = map { $_ => 1 } @robin;
            foreach my $dir (sort keys %pool_dirs) {
                push(@robin, $dir) unless exists $current{$dir};
            }
            my $max_robin = @robin;
            
            my $req = $self->new_requirements(memory => 100, time => 1);
            foreach my $file (@{$self->inputs->{files}}) {
                # pick a pool directory based on round-robin
                my $root_dir;
                do {
                    $root_dir = $robin[$robin_index++];
                    $robin_index %= $max_robin;
                } while (! exists $dirs{$root_dir});
                
                # pick subdirs based on md5sum of $file's id
                my $dmd5 = Digest::MD5->new();
                my $ifile_id = $file->id;
                $dmd5->add($ifile_id);
                my $md5 = $dmd5->hexdigest;
                my @chars = split("", $md5);
                
                my $ofile = $self->output_file(output_key => 'moved_files',
                                               basename => $file->basename,
                                               output_dir => dir($root_dir, @chars[0..3], $ifile_id),
                                               type => $file->type,
                                               metadata => $file->metadata)->path;
                
                # VRPipe::File->move will do a copy, check md5s, and only then
                # delete the source, so it is nice and safe
                my $ifile = $file->path;
                $self->dispatch_vrpipecode("VRPipe::File->get(path => q[$ifile])->move(VRPipe::File->get(path => q[$ofile]), check_md5s => 1);", $req);
            }
        };
    }
    method outputs_definition {
        return { moved_files => VRPipe::StepIODefinition->get(type => 'any',
                                                              description => 'files at their new locations',
                                                              max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Safely move files from one disc to a pool of one or more other discs, eg. for archival purposes";
    }
    method max_simultaneous {
        return 25;
    }
}

1;
