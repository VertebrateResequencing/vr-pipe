=head1 NAME

VRPipe::Steps::bwa_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::bwa_index with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file to map against'),
                 bwa_index_options => VRPipe::StepOption->get(description => 'options to bwa index, excluding the reference fasta file',
                                                              optional => 1,
                                                              default_value => '-a bwtsw'),
                 bwa_exe => VRPipe::StepOption->get(description => 'path to your bwa executable',
                                                    optional => 1,
                                                    default_value => 'bwa') };
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
            
            my $bwa_exe = $options->{bwa_exe};
            my $bwa_opts = $options->{bwa_index_options};
            if ($bwa_opts =~ /$ref|index/) {
                $self->throw("bwa_index_options should not include the reference or index subcommand");
            }
            my $cmd = $bwa_exe.' index '.$bwa_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->get(exe => 'bwa', version => VRPipe::StepCmdSummary->determine_version($bwa_exe, '^Version: (.+)$'), summary => 'bwa index '.$bwa_opts.' $reference_fasta'));
            $cmd .= ' '.$ref;
            
            foreach my $suffix (qw(bwt pac rbwt rpac rsa sa)) {
                $self->output_file(output_key => 'bwa_index_binary_files', output_dir => $ref->dir->stringify, basename => $ref->basename.'.'.$suffix, type => 'bin');
            }
            foreach my $suffix (qw(amb ann)) {
                $self->output_file(output_key => 'bwa_index_text_files', output_dir => $ref->dir->stringify, basename => $ref->basename.'.'.$suffix, type => 'txt');
            }
            
            $self->dispatch([$cmd, $self->new_requirements(memory => 3900, time => 1), {block_and_skip_if_ok => 1}]);
        };
    }
    method outputs_definition {
        return { bwa_index_binary_files => VRPipe::StepIODefinition->get(type => 'bin', description => 'the files produced by bwa index', min_files => 6, max_files => 6),
                 bwa_index_text_files => VRPipe::StepIODefinition->get(type => 'txt', description => 'the files produced by bwa index', min_files => 2, max_files => 2) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Indexes a reference genome fasta file, making it suitable for use in subsequent bwa mapping";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
