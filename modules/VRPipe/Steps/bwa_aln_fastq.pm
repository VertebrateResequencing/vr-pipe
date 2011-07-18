use VRPipe::Base;

class VRPipe::Steps::bwa_aln_fastq with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file to map against'),
                 bwa_aln_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa aln command line, including desired options, but excluding the input fastq, reference and -f option',
                                                        optional => 1,
                                                        default_value => 'bwa aln -q 15') };
    }
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->get(type => 'fq',
                                                              max_files => -1,
                                                              description => 'fastq files, which will be alnd independently',
                                                              metadata => {reads => 'total number of reads (sequences)',
                                                                           paired => '0=unpaired; 1=reads in this file are forward; 2=reads in this file are reverse',
                                                                           mate => 'if paired, the path to the fastq that is our mate',
                                                                           optional => ['mate']}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $ref = Path::Class::File->new($self->options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $cmd = $self->options->{bwa_aln_cmd};
            my @c = split(" ", $cmd);
            unless ($c[1] eq 'aln') {
                $self->throw("bad bwa_aln_cmd '$cmd'");
            }
            
            if ($cmd =~ /$ref|-f/) {
                $self->throw("bwa_aln_cmd should not include the reference or -f option");
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $fastq (@{$self->inputs->{fastq_files}}) {
                my $sai_file = $self->output_file(output_key => 'bwa_sai_files',
                                                  output_dir => $fastq->dir,
                                                  basename => $fastq->basename.'.sai',
                                                  type => 'bin',
                                                  metadata => {source_fastq => $fastq->path->stringify,
                                                               reference_fasta => $ref->stringify});
                my $this_cmd = $cmd.' -f '.$sai_file->path.' '.$ref.' '.$fastq->path;
                $self->dispatch([$this_cmd, $req]);
            }
        };
    }
    method outputs_definition {
        return { bwa_sai_files => VRPipe::StepIODefinition->get(type => 'bin',
                                                                max_files => -1,
                                                                description => 'output files of independent bwa aln calls on each input fastq',
                                                                metadata => {source_fastq => 'the fastq file that was input to bwa aln to generate this sai file',
                                                                             reference_fasta => 'the reference fasta that reads were aligned to'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Aligns the input fastq(s) with bwa to the reference";
    }
}

1;
