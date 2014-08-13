
=head1 NAME

VRPipe::Steps::samtools_fasta_gc_stats - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::samtools_fasta_gc_stats extends VRPipe::Steps::plot_bamstats {
    around options_definition {
        return {
            %{ $self->$orig },
            reference_fasta    => VRPipe::StepOption->create(description => 'absolute path to genome reference file'),
            exome_targets_file => VRPipe::StepOption->create(
                description => 'absolute path to a file describing the targets/baits used for exome pulldown (tab-delimited [chr,start,end], where start is 1-based, and end is inclusive)',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $options     = $self->options;
            my $plot_bs_exe = $options->{plot_bamstats_exe};
            
            # we don't accept opts from user, since we handle them all ourselves
            # and they must not be overridden
            my $ref = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            my $plot_bs_opts = '-s ' . $ref;
            
            my $output_basename = $ref->basename . '.gc_stats';
            
            my $targets = $options->{exome_targets_file};
            if ($targets) {
                $targets = file($targets);
                $self->throw("exome_targets_file must be an absolute path") unless $targets->is_absolute;
                $plot_bs_opts .= ' -t ' . $targets;
                
                # we need to uniqify the basename for each target file
                my $tfile = VRPipe::File->create(path => $targets);
                unless ($tfile->md5) {
                    $tfile->update_md5;
                }
                $output_basename .= '.targeted-' . $tfile->md5;
            }
            
            my $output_path = $self->output_file(output_key => 'fasta_gc_stats_file', output_dir => $ref->dir, basename => $output_basename, type => 'txt')->path;
            
            my $cmd = qq[$plot_bs_exe $plot_bs_opts > $output_path];
            $self->dispatch([$cmd, $self->new_requirements(memory => 500, time => 1), { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return {
            fasta_gc_stats_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that describes the gc stats of an input fasta'
            )
        };
    }
    
    method description {
        return "Uses plot-bamcheck to generate a file of gc stats for a fasta file";
    }
}

1;
