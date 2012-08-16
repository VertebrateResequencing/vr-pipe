
=head1 NAME

VRPipe::Steps::plot_bamcheck - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Steps::plot_bamcheck with VRPipe::StepRole {
    method options_definition {
        return {
            plot_bamcheck_exe => VRPipe::StepOption->create(
                description   => 'path to plot-bamcheck executable',
                optional      => 1,
                default_value => 'plot-bamcheck'
            )
        };
    }
    
    method inputs_definition {
        return {
            bamcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'bamcheck output files',
                max_files   => -1,
                metadata    => {
                    source_bam => 'path to the bam file this bamcheck file was created from',
                    lane       => 'lane name (a unique identifer for this sequencing run, aka read group)'
                }
            ),
            fasta_gc_stats_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'plot-bamcheck -s output file for the reference fasta')
        };
    }
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $options     = $self->options;
            my $req         = $self->new_requirements(memory => 500, time => 1);
            my $plot_bc_exe = $options->{plot_bamcheck_exe};
            
            # we don't accept opts from user, since we handle them all ourselves
            # and they must not be overridden
            my $plot_bc_opts = '-r ' . $self->inputs->{fasta_gc_stats_file}->[0]->path;
            
            foreach my $bc_file (@{ $self->inputs->{bamcheck_files} }) {
                my $meta       = $bc_file->metadata;
                my $prefix     = $meta->{lane} || $bc_file->basename;
                my $source_bam = $meta->{source_bam};
                my $sb_meta    = VRPipe::File->get(path => $source_bam)->metadata;
                
                $self->output_file(temporary  => 1,                basename => $prefix . '-quals.gp',     type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-quals.png',    type => 'bin', metadata => { caption => 'Qualities', source_bam => $source_bam });
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-quals2.png',   type => 'bin', metadata => { caption => 'Qualities', source_bam => $source_bam });
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-quals3.png',   type => 'bin', metadata => { caption => 'Qualities', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-quals-hm.gp',  type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-quals-hm.png', type => 'bin', metadata => { caption => 'Qualities', source_bam => $source_bam });
                if ((exists $sb_meta->{mean_insert_size} && $sb_meta->{mean_insert_size} > 0) || (exists $sb_meta->{targeted_mean_insert_size} && $sb_meta->{targeted_mean_insert_size} > 0)) {
                    $self->output_file(temporary => 1, basename => $prefix . '-insert-size.gp', type => 'txt');
                    $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-insert-size.png', type => 'bin', metadata => { caption => 'Insert Size', source_bam => $source_bam });
                }
                $self->output_file(temporary  => 1,                basename => $prefix . '-gc-content.gp',      type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-gc-content.png',     type => 'bin', metadata => { caption => 'GC Content', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-gc-depth.gp',        type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-gc-depth.png',       type => 'bin', metadata => { caption => 'GC Depth', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-acgt-cycles.gp',     type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-acgt-cycles.png',    type => 'bin', metadata => { caption => 'ACGT Cycles', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-coverage.gp',        type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-coverage.png',       type => 'bin', metadata => { caption => 'Coverage', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-indel-dist.gp',      type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-indel-dist.png',     type => 'bin', metadata => { caption => 'Indel distribution', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-indel-cycles.gp',    type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-indel-cycles.png',   type => 'bin', metadata => { caption => 'Indels per cycle', source_bam => $source_bam });
                $self->output_file(temporary  => 1,                basename => $prefix . '-mism-per-cycle.gp',  type => 'txt');
                $self->output_file(output_key => 'bamcheck_plots', basename => $prefix . '-mism-per-cycle.png', type => 'bin', metadata => { caption => 'Mismatches per cycle', source_bam => $source_bam });
                
                my $bc_path = $bc_file->path;
                my $cmd     = qq[$plot_bc_exe $plot_bc_opts -p $prefix $bc_path];
                $self->dispatch([$cmd, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bamcheck_plots => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'png files produced by plot-bamcheck, with a caption in the metadata',
                min_files   => 11,
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Uses plot-bamcheck to draw plots from the associated bamcheck file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
