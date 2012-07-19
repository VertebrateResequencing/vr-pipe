use VRPipe::Base;

class VRPipe::Steps::test_step_fail with VRPipe::StepRole {
    method options_definition {
        return {};
    }
    
    method inputs_definition {
        return { file => VRPipe::StepIODefinition->create(type => 'any', description => 'an in file') };
    }
    
    method body_sub {
        return sub {
            my $self   = shift;
            my $req    = $self->new_requirements(memory => 50, time => 1);
            my ($file) = @{ $self->inputs->{file} };
            my $out    = $self->output_file(basename => $file->basename, output_key => 'out', type => 'txt')->path;
            $self->dispatch_vrpipecode("use VRPipe::Steps::test_step_fail; VRPipe::Steps::test_step_fail->failer(out => q[$out]);", $req);
        };
    }
    
    method outputs_definition {
        return { out => VRPipe::StepIODefinition->create(type => 'txt', description => 'an out file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "test step that fails twice before working";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method failer (ClassName|Object $self: Str :$out!) {
        my $made;
        my @progress_files;
        foreach my $i (1 .. 3) {
            my $ofile = VRPipe::File->create(path => $out . '.' . $i);
            push(@progress_files, $ofile);
            next if $ofile->s;
            my $ofh = $ofile->openw;
            print $ofh "$i\n";
            $ofile->close;
            $made = $i;
            last;
        }
        
        if ($made == 3) {
            my $ofile = VRPipe::File->get(path => $out);
            my $ofh = $ofile->openw;
            print $ofh "success\n";
            $ofile->close;
            foreach my $pfile (@progress_files) {
                $pfile->unlink;
            }
        }
        else {
            print "stdout message: failing on purpose since this is try $made\n";
            die "stderr message: failing on purpose since this is try $made\n";
        }
    }
}

1;
