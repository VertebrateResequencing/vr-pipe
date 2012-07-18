
=head1 NAME

VRPipe::Steps::stampy_buildhash - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::stampy_buildhash with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file to map against'),
                 stampy_exe      => VRPipe::StepOption->create(
                                                          description   => 'path to your stampy.py executable',
                                                          optional      => 1,
                                                          default_value => 'stampy.py') };
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
            
            my $stampy_exe = $options->{stampy_exe};
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'stampy', version => VRPipe::StepCmdSummary->determine_version($stampy_exe, '^stampy v(\S+)'), summary => 'stampy.py -g $ref.fa -H $ref.fa'));
            
            my $cmd = $stampy_exe . " -g $ref -H $ref";
            $self->output_file(output_key => 'stampy_index_sthash_file',
                               output_dir => $ref->dir->stringify,
                               basename   => $ref->basename . '.sthash',
                               type       => 'bin');
            $self->dispatch([$cmd, $self->new_requirements(memory => 3900, time => 1), { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return { stampy_index_sthash_file => VRPipe::StepIODefinition->create(type => 'bin', description => 'the file produced by stampy -H') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Part 2 (of 2) of indexing a reference genome fasta file, making it suitable for use in subsequent stampy mapping";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
