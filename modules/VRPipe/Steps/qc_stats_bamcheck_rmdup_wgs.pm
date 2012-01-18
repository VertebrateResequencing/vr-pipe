use VRPipe::Base;
use Data::Dumper;
class VRPipe::Steps::qc_stats_bamcheck_rmdup_wgs with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->get(description => 'path to your bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'bamcheck -d'),
               };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', description => 'bam files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $bamcheck_exe = $options->{bamcheck_exe};
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $ifile = $bam_file->path;
                $self->output_file(output_key => 'bam_files_with_metadata', output_dir => $ifile->dir, basename => $ifile->basename, type => 'bam');
                my $meta = $bam_file->metadata; 
                unless ($meta->{rmdup_bases} && $meta->{rmdup_reads}) {
                    my $check_file = $self->output_file(output_key => 'bamcheck_files', basename => $ifile->basename.'.bamcheck', type => 'txt');
                    my $ofile = $check_file->path;
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::qc_stats_bamcheck_rmdup_wgs', 'stats_from_bamcheck', ["$bamcheck_exe $ifile > $ofile", $req, {output_files => [$check_file]}]);
                }
            }
        };
    }

    method outputs_definition {
         return { bamcheck_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                                 description => 'bamcheck file required for graph step',
                                                                 max_files => -1,
                                                                 metadata => { sample => 'sample name' }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Produces a bamcheck output file that has metadata added to it from the bam file that is used to draw plots in the next step";
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
            my $check_meta = {};
            # parse the bamcheck file
            my $parser = VRPipe::Parser->create('bamcheck', {file => $check_file});
            open(my $ifh, $check_path) || die "could not open $check_path\n";
    		$new_meta->{rmdup_reads} = $parser->sequences;
			$new_meta->{rmdup_reads_mapped} = $parser->reads_mapped;
			$new_meta->{rmdup_bases_mapped_c} = $parser->bases_mapped_cigar;
			$new_meta->{rmdup_bases} = $parser->total_length;
			$new_meta->{rmdup_bases_trimmed} = $parser->bases_trimmed;
            $bam_file->add_metadata($new_meta);
            my $meta = $bam_file->metadata;
            $check_meta->{sample} = $meta->{sample};
            $check_file->add_metadata($check_meta);
        }
        else {
           $self->throw("$check_path failed to be made");
        }
    }
}

1;
