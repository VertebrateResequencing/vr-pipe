
=head1 NAME

VRPipe::Steps::smalt_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::smalt_index with VRPipe::StepRole {
    method options_definition {
        return {
            reference_fasta     => VRPipe::StepOption->create(description => 'absolute path to genome reference file to map against'),
            smalt_index_options => VRPipe::StepOption->create(
                description   => 'options to bwa index, excluding the reference fasta file',
                optional      => 1,
                default_value => '-k 13 -s 4'
            ),
            smalt_exe => VRPipe::StepOption->create(description => 'path to your smalt executable', optional => 1, default_value => 'smalt')
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
            
            my $smalt_exe  = $options->{smalt_exe};
            my $smalt_opts = $options->{smalt_index_options};
            if ($smalt_opts =~ /$ref|index/) {
                $self->throw("smalt_index_options should not include the reference or index subcommand");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'smalt',
                    version => VRPipe::StepCmdSummary->determine_version($smalt_exe . ' version', '^Version: (.+)$'),
                    summary => 'smalt index ' . $smalt_opts . ' $index_base $reference_fasta'
                )
            );
            
            my $index_base = $smalt_opts;
            $index_base =~ s/[^ks0-9]//g;
            $index_base = $ref->basename . '_' . $index_base;
            my $index_base_path = Path::Class::File->new($ref->dir, $index_base);
            my $cmd = $smalt_exe . ' index ' . $smalt_opts . ' ' . $index_base_path . ' ' . $ref;
            
            foreach my $suffix (qw(sma smi)) {
                $self->output_file(
                    output_key => 'smalt_index_binary_files',
                    output_dir => $ref->dir->stringify,
                    basename   => $index_base . '.' . $suffix,
                    type       => 'bin',
                    metadata   => { index_base => $index_base_path->stringify }
                );
            }
            
            $self->dispatch([$cmd, $self->new_requirements(memory => 6000, time => 1), { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return {
            smalt_index_binary_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'the files produced by smalt index',
                min_files   => 2,
                max_files   => 2,
                metadata    => { index_base => 'base used in smalt index command' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Indexes a reference genome fasta file, making it suitable for use in subsequent smalt mapping";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
