use VRPipe::Base;

class VRPipe::Steps::bwa_aln_fastq with VRPipe::StepRole {
    method options_definition {
        return { bwa_aln_cmd => VRPipe::StepOption->get(description => 'the near-complete bwa aln command line, including desired options, but excluding the input fastq',
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
            $self->throw("foo");
            #my ($lane, $se, @pe, %ended);
            #foreach my $fastq (@{$self->inputs->{fastq_files}}) {
            #    my $metadata = $fastq->metadata;
            #    $lane ||= $metadata->{lane};
            #    my $paired = $metadata->{paired};
            #    if ($paired == 0) {
            #        $se = $fastq;
            #    }
            #    elsif ($paired == 1) {
            #        unshift(@pe, $fastq);
            #    }
            #    elsif ($paired == 2) {
            #        push(@pe, $fastq);
            #    }
            #}
            #$ended{pe} = @pe == 2 ? \@pe : undef;
            #$ended{se} = $se ? [$se] : undef;
            #my $out_root = dir($self->output_root, $lane);
        };
    }
    method outputs_definition {
        return { bwa_sai_files => VRPipe::StepIODefinition->get(type => 'dat',
                                                                max_files => -1,
                                                                description => 'output files of independent bwa aln calls on each input fastq',
                                                                metadata => {source_fastq => 'the fastq file that was input to bwa aln to generate this sai file'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Aligns the input fastq(s) with bwa to the reference";
    }
}

1;
