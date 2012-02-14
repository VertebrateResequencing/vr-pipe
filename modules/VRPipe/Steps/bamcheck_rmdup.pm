use VRPipe::Base;

class VRPipe::Steps::bamcheck_rmdup extends VRPipe::Steps::bamcheck {
    use VRPipe::Parser;
        method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->get(description => 'path to your bamcheck executable',
                                                         optional => 1,
                                                         default_value => 'bamcheck'),
                 bamcheck_rmdup_options => VRPipe::StepOption->get(description => 'options to bamcheck, excluding the value to -r which will come from reference_fasta option',
                                                         optional => 1,
                                                         default_value => '-d')};
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $bamcheck_exe = $options->{bamcheck_exe};
            my $bc_rmdup_opts = $options->{bamcheck_rmdup_options};
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $ifile = $bam_file->path;
                $self->output_file(output_key => 'bam_files_with_metadata', output_dir => $ifile->dir, basename => $ifile->basename, type => 'bam');
                my $check_file = $self->output_file(basename => $ifile->basename.'.rmdup.bamcheck', type => 'txt', temporary => 1);
                my $ofile = $check_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bamcheck', 'stats_from_bamcheck', ["$bamcheck_exe $bc_rmdup_opts $ifile > $ofile", $req, {output_files => [$check_file]}]);
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
                                                                                       optional => ['library', 'sample', 'center_name', 'platform', 'study', 'insert_size']}),
                   bamcheck_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'bamcheck output files', max_files => -1) 
               }; 
    }
    method description {
        return "Creates a bamcheck file with duplicates removed and associates some of the stats as metadata on the bam file";
    }

}

1;
