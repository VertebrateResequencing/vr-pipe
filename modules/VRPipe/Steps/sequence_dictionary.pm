
=head1 NAME

VRPipe::Steps::sequence_dictionary - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Steps::sequence_dictionary with VRPipe::StepRole {
    use Digest::MD5 qw(md5_hex);
    use VRPipe::Parser;
    
    method _build_smaller_recommended_requirements_override {
        return 0;
    }
    
    method options_definition {
        return {
            reference_fasta         => VRPipe::StepOption->create(description => 'absolute path to fasta file'),
            reference_assembly_name => VRPipe::StepOption->create(description => 'public name of the assembly, eg. NCBI37; defaults to being excluded', optional => 1),
            reference_public_url    => VRPipe::StepOption->create(description => 'public url that the reference_fasta can be accessed from; defaults to reference_fasta path', optional => 1),
            reference_species       => VRPipe::StepOption->create(description => 'species of the reference genome; defaults to being excluded', optional => 1)
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $ref     = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $dict_file = $self->output_file(output_key => 'reference_dict', output_dir => $ref->dir->stringify, basename => $ref->basename . '.dict', type => 'txt')->path;
            my $ur = $options->{reference_public_url} || 'file:' . $ref;
            my $as = $options->{reference_assembly_name};
            my $sp = $options->{reference_species};
            
            my @constants = ('UR:' . $ur);
            push(@constants, 'AS:' . $as) if $as;
            push(@constants, 'SP:' . $sp) if $sp;
            my $constants = join("\t", @constants);
            
            my $code = qq[use VRPipe::Steps::sequence_dictionary; VRPipe::Steps::sequence_dictionary->dicter(ref => q[$ref], dict => q[$dict_file], constants => q[$constants]);];
            
            $self->dispatch_vrpipecode($code, $self->new_requirements(memory => 2000, time => 1), { block_and_skip_if_ok => 1 });
        };
    }
    
    method outputs_definition {
        return { reference_dict => VRPipe::StepIODefinition->create(type => 'txt', description => '.dict file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Creates a sequence dictionary (.dict file) from a fasta file, suitable for use by Picard/GATK and for forming good bam headers";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method dicter (ClassName|Object $self: Str|File :$ref!, Str|File :$dict!, Str :$constants?) {
        my $pars = VRPipe::Parser->create('fasta', { file => $ref });
        
        my $dict_file = VRPipe::File->get(path => $dict);
        $dict_file->disconnect;
        
        my $dict_content = "\@HD\tVN:1.0\tSO:unsorted\n";
        
        my $pr = $pars->parsed_record();
        while ($pars->next_record()) {
            my $sn = 'SN:' . $pr->[0];
            
            my $seq = uc($pr->[1]);
            $seq =~ s/\s//g;
            my $ln = 'LN:' . length($seq);
            my $m5 = 'M5:' . md5_hex($seq);
            
            my @fields = ('@SQ', $sn, $ln, $m5);
            push(@fields, $constants) if ($constants && $constants =~ /\S/);
            
            $dict_content .= join("\t", @fields) . "\n";
        }
        
        my $write = 1;
        if (-s $dict_file->path) {
            # There is an obscure issue possibly related to race conditions
            # where the dict file can have size already, even though it never
            # existing before this pipelinesetup wanted to create it.
            # Presumably another submission for this same job called a run
            # a few miliseconds before us. We don't want to throw here, however,
            # or we'll end up with a complete .dict file but a job status of
            # failed, then subsequent retries will all fail.
            # Instead we'll compare to what the existing .dict file contains and
            # fail/succeed accordingly
            sleep(5);
            my $ifh             = $dict_file->openr;
            my $current_content = $dict_file->slurp;
            unless ($current_content eq $dict_content) {
                $self->warn("A dict file with different content already exists at '" . $dict_file->path . "'; was it made by a different pipeline with different settings? I will not overwrite it");
                $write = 0;
            }
        }
        
        if ($write) {
            my $ofh = $dict_file->openw;
            print $ofh $dict_content;
            $dict_file->close;
        }
        
        # GATK adds a dict file automatically when it runs if it doesn't already exist
        # Unformunately it replaces .fa or .fasta with .dict rather than appending.
        # Create a symlink to our just created file to replace this. This will prevent
        # runs where many GATK jobs start at once and all want to create this file.
        my $symlink_path = $ref;
        $symlink_path =~ s/fa(sta)?(\.gz)?$/dict/;
        my $dict_symlink = VRPipe::File->create(path => $symlink_path);
        if (-s $symlink_path && !(-l $symlink_path)) {
            $dict_symlink->unlink;
        }
        $dict_file->symlink($dict_symlink);
    }
}

1;
