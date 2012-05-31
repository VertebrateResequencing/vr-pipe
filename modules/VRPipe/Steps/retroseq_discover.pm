use VRPipe::Base;

class VRPipe::Steps::retroseq_discover with VRPipe::StepRole {
    method options_definition {
        return { 
            retroseq_exe => VRPipe::StepOption->get(description => 'full path to retroseq.pl', optional => 1, default_value => 'retroseq.pl'),
            refTEs_param => VRPipe::StepOption->get(description => '-refTEs option value, tab file with TE type and BED file of reference elements', optional => 1),
            exde_param => VRPipe::StepOption->get(description => '-exd option value, fofn of BED files of discordant mate exclusion regions', optional => 1),
        };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                            description => 'bam files',
                                                            max_files => -1) };
    }

    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            my $retroseq_exe = $options->{retroseq_exe};
            my $refTEs_param = $options->{refTEs_param};
            my $exde_param = $options->{exde_param};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $basename = $bam_file->basename;
                $basename =~ s/\.bam$/.cand.tab/;
                my $rseq_bed = $self->output_file(output_key => 'rseq_bed', basename => $basename, type => 'txt',
													metadata => {source_bam => $bam_file->path->stringify});
                
                my $input_path = $bam_file->path;
                my $output_path = $rseq_bed->path;
                
                my $cmd = "$retroseq_exe -discover -bam $input_path";
                $cmd .= " -refTEs $refTEs_param" if $refTEs_param;
                $cmd .= " -exd $exde_param" if $exde_param;
                $cmd .= " -output $output_path";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::retroseq_discover', 'gen_candidates', [$cmd, $req, {output_files => [$rseq_bed]}]);
            }
        };
    }
    method outputs_definition {
        return { rseq_bed => VRPipe::StepIODefinition->get(type => 'txt',
                                                               description => 'retroseq -discover output in BED format',
                                                               max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generates retroseq candidate supporting read pairs in BED format from input bam";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method gen_candidates (ClassName|Object $self: Str $cmd_line) {

        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ /-output (\S+)$/;
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        
        if ($output_file->lines == 0) {
            $output_file->unlink;
			$self->throw("Output $output_path is empty)");
        }
        else {
            return 1;
        }
    }
}

1;
