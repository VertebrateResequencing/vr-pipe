use VRPipe::Base;

class VRPipe::Steps::bam_reheader with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->get(description => 'path to your samtools executable',
                                                         optional => 1,
                                                         default_value => 'samtools') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                            max_files => -1,
                                                            description => '1 or more bam files',
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                         library => 'library name',
                                                                         sample => 'sample name',
                                                                         center_name => 'center name',
                                                                         platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                         study => 'name of the study',
                                                                         insert_size => 'expected (mean) insert size if paired',,
                                                                         bases => 'total number of base pairs',
                                                                         reads => 'total number of reads (sequences)',
                                                                         paired => '0=unpaired reads were mapped; 1=paired reads were mapped',
                                                                         optional => ['library', 'insert_size', 'sample', 'center_name', 'platform', 'study']}),
                 dict_file => VRPipe::StepIODefinition->get(type => 'txt',
                                                            description => 'a sequence dictionary file for your reference fasta') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            my $dict_path = $self->inputs->{dict_file}->[0]->path;
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            my $step_state = $self->step_state->id;
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam->path;
                my $bam_meta = $bam->metadata;
                my $basename = $bam->basename;
                my $headed_bam_file = $self->output_file(output_key => 'headed_bam_files',
                                                  basename => $basename,
                                                  type => 'bam',
                                                  metadata => $bam_meta);
                
                my $headed_bam_path = $headed_bam_file->path;
                
                $self->output_file(output_key => 'temp_header_file',
                                     basename => $basename.'.header',
                                     type => 'txt',
                                     temporary => 1);

                my $this_cmd = "use VRPipe::Steps::bam_reheader; VRPipe::Steps::bam_reheader->reheader_and_check(samtools => q[$samtools], dict => q[$dict_path], output => q[$headed_bam_path], step_state => $step_state, bam => q[$bam_path]);";
                $self->dispatch_vrpipecode($this_cmd, $req); # deliberately do not include {output_files => [$headed_bam_file]} so that any temp files we made will get their stats updated prior to auto-deletion
            }
        };
    }
    method outputs_definition {
        return { headed_bam_files => VRPipe::StepIODefinition->get(type => 'bam',
                                                          max_files => -1,
                                                          description => 'a bam file with good header',
                                                          metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)',
                                                                       library => 'library name',
                                                                       sample => 'sample name',
                                                                       center_name => 'center name',
                                                                       platform => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                                                                       study => 'name of the study, put in the DS field of the RG header line',
                                                                       insert_size => 'expected (mean) insert size if paired',,
                                                                       bases => 'total number of base pairs',
                                                                       reads => 'total number of reads (sequences)',
                                                                       paired => '0=unpaired reads were mapped; 1=paired reads were mapped; 2=mixture of paired and unpaired reads were mapped',
                                                                       optional => ['library', 'insert_size', 'sample', 'center_name', 'platform', 'study']}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Replaces a bam header so that it has complete sequence information, a good RG line, and chained PG lines";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    method reheader_and_check (ClassName|Object $self: Str|File :$samtools, Str|File :$dict, Str|File :$output, Persistent :$step_state, Str|File :$bam) {
        # make a nice sam header
        my $header_file = VRPipe::File->get(path => $output.'.header');
        my $header_path = $header_file->path;
        my $hfh = $header_file->openw;
        
        my $header_lines = 0;
        print $hfh "\@HD\tVN:1.0\tSO:coordinate\n";
        $header_lines++;
        
        # copy over the SQ lines from the dict file
        my $dict_file = VRPipe::File->get(path => $dict);
        my $dfh = $dict_file->openr;
        while (<$dfh>) {
            next unless /^\@SQ/;
            print $hfh $_;
            $header_lines++;
        }
        
        # construct the RG line from the bam metadata
        my $headed_bam_file = VRPipe::File->get(path => $output);
        my $meta = $headed_bam_file->metadata;
        print $hfh "\@RG\tID:", $meta->{lane};
        if (defined $meta->{library}) {
            print $hfh "\tLB:", $meta->{library};
        }
        if (defined $meta->{sample}) {
            print $hfh "\tSM:", $meta->{sample};
        }
        if (defined $meta->{insert_size}) {
            print $hfh "\tPI:", $meta->{insert_size};
        }
        if (defined $meta->{center_name}) {
            print $hfh "\tCN:", $meta->{center_name};
        }
        if (defined $meta->{platform}) {
            print $hfh "\tPL:", $meta->{platform};
        }
        if (defined $meta->{study}) {
            print $hfh "\tDS:", $meta->{study};
        }
        print $hfh "\n";
        $header_lines++;
        
        # construct a chain of PG lines for the header by looking at previous
        # steps in our pipeline
        my $this_step_state = VRPipe::StepState->get(id => $step_state);
        my $pipelinesetup = $this_step_state->pipelinesetup;
        my $dataelement = $this_step_state->dataelement;
        my $stepmember = $this_step_state->stepmember;
        my $this_stepm_id = $stepmember->id;
        my $pipeline = $stepmember->pipeline;
        my $pp;
        foreach my $stepm ($pipeline->steps) {
            last if $stepm->id == $this_stepm_id;
            
            my $cmd_summary = VRPipe::StepState->get(pipelinesetup => $pipelinesetup, stepmember => $stepm, dataelement => $dataelement)->cmd_summary || next;
            
            my $step_name = $stepm->step->name;
            print $hfh "\@PG\tID:$step_name\tPN:", $cmd_summary->exe, "\t";
            if ($pp) {
                print $hfh "PP:$pp\t";
            }
            print $hfh "VN:", $cmd_summary->version, "\tCL:", $cmd_summary->summary, "\n";
            $header_lines++;
            
            $pp = $step_name;
        }
        $header_file->close;
        
        my $cmd_line = qq[$samtools reheader $header_path $bam > $output];
        $headed_bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $expected_lines = $meta->{reads} + $header_lines;
        $headed_bam_file->update_stats_from_disc(retries => 3);
        my $actual_lines = $headed_bam_file->lines;
        
        if ($actual_lines == $expected_lines) {
            return 1;
        }
        else {
            $headed_bam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_lines lines were generated in the reheaded bam file, yet there were $expected_lines records in the input bam files and header");
        }
    }
}

1;
