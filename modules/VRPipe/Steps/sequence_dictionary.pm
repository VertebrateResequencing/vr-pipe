use VRPipe::Base;

class VRPipe::Steps::sequence_dictionary with VRPipe::StepRole {
    use Digest::MD5 qw(md5_hex);
    use VRPipe::Parser;
    
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to fasta file'),
                 reference_assembly_name => VRPipe::StepOption->get(description => 'public name of the assembly, eg. NCBI37; defaults to being excluded', optional => 1),
                 reference_public_url => VRPipe::StepOption->get(description => 'public url that the reference_fasta can be accessed from; defaults to reference_fasta path', optional => 1),
                 reference_species => VRPipe::StepOption->get(description => 'species of the reference genome; defaults to being excluded', optional => 1)};
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $dict_file = $self->output_file(output_key => 'reference_dict', output_dir => $ref->dir->stringify, basename => $ref->basename.'.dict', type => 'txt')->path;
            my $ur = $options->{reference_public_url} || 'file:'.$ref;
            my $as = $options->{reference_assembly_name};
            my $sp = $options->{reference_species};
            
            my @constants = ('UR:'.$ur);
            push(@constants, 'AS:'.$as) if $as;
            push(@constants, 'SP:'.$sp) if $sp;
            my $constants = join("\t", @constants);
            
            my $code = qq[use VRPipe::Steps::sequence_dictionary; VRPipe::Steps::sequence_dictionary->dicter(ref => q[$ref], dict => q[$dict_file], constants => q[$constants]);];
            
            $self->dispatch_vrpipecode($code, $self->new_requirements(memory => 2000, time => 1), {block_and_skip_if_ok => 1});
        };
    }
    method outputs_definition {
        return { reference_dict => VRPipe::StepIODefinition->get(type => 'txt', description => '.dict file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Creates a sequence dictionary (.dict file) from a fasta file, suitable for use by Picard/GATK and for forming good bam headers";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method dicter (ClassName|Object $self: Str|File :$ref!, Str|File :$dict!, Str :$constants?) {
        my $pars = VRPipe::Parser->create('fasta', {file => $ref});
        
        my $dict_file = VRPipe::File->get(path => $dict);
        my $ofh = $dict_file->openw;
        $dict_file->disconnect;
        
        print $ofh "\@HD\tVN:1.0\tSO:unsorted\n";
        
        my $pr = $pars->parsed_record();
        while ($pars->next_record()) {
            my $sn = 'SN:'.$pr->[0];
            
            my $seq = uc($pr->[1]);
            $seq =~ s/\s//g;
            my $ln = 'LN:'.length($seq);
            my $m5 = 'M5:'.md5_hex($seq);
            
            my @fields = ('@SQ', $sn, $ln, $m5);
            push(@fields, $constants) if ($constants && $constants =~ /\S/);
            
            print $ofh join("\t", @fields), "\n";
        }
        
        $dict_file->close;
    }
}

1;
