
=head1 NAME

VRPipe::Steps::sga_generate_sampled_suffix_array - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

# Usage: sga gen-ssa [OPTION] ... READSFILE
# Build a sampled suffix array for the reads in READSFILE using the BWT
#
#   -v, --verbose                        display verbose output
#       --help                           display this help and exit
#   -t, --threads=NUM                    use NUM threads to construct the index (default: 1)
#   -c, --check                          validate that the suffix array/bwt is correct

use VRPipe::Base;

class VRPipe::Steps::sga_generate_sampled_suffix_array with VRPipe::StepRole {
    method options_definition {
        return {
            sga_gen_ssa_options => VRPipe::StepOption->create(description => 'options to sga gen-ssa',      optional => 1),
            sga_exe             => VRPipe::StepOption->create(description => 'path to your sga executable', optional => 1, default_value => 'sga')
        };
    }
    
    method inputs_definition {
        return { reference_fasta => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'reference fasta file') },;
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $ref_file = $self->inputs->{reference_fasta}[0];
            
            my $sga_exe  = $options->{sga_exe};
            my $sga_opts = $options->{sga_gen_ssa_options};
            if ($sga_opts =~ /gen-ssa/) {
                $self->throw("sga_gen_ssa_options should not include the reference or gen-ssa subcommand");
            }
            my $cmd = $sga_exe . ' gen-ssa ' . $sga_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga gen-ssa ' . $sga_opts . ' $reference_fasta'));
            $cmd .= ' ' . $ref_file->path;
            
            my $basename = $ref_file->basename;
            $basename =~ s/(fa|fasta)(\.gz)?$/ssa/;
            my $ssa_file = $self->output_file(output_key => 'sampled_suffix_array', output_dir => $ref_file->dir->stringify, basename => $basename, type => 'bin');
            
            my $req = $self->new_requirements(memory => 9000, time => 1);
            $self->dispatch([$cmd, $req, { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return { sampled_suffix_array => VRPipe::StepIODefinition->create(type => 'bin', description => 'the sampled suffix array file', max_files => 1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Build a sampled suffix array for the reads in READSFILE using the BWT";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
