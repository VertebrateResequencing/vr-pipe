
=head1 NAME

VRPipe::Steps::plot_bamstats - a step

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

class VRPipe::Steps::plot_bamstats with VRPipe::StepRole {
    use VRPipe::Schema;
    
    method options_definition {
        return {
            plot_bamstats_exe => VRPipe::StepOption->create(
                description   => 'path to plot-bamstats executable',
                optional      => 1,
                default_value => 'plot-bamstats'
            )
        };
    }
    
    method inputs_definition {
        return {
            stats_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'samtools stats output files',
                max_files   => -1
            ),
            fasta_gc_stats_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'plot-bamstats -s output file for the reference fasta')
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $req           = $self->new_requirements(memory => 500, time => 1);
            my $plot_bamstats = $options->{plot_bamstats_exe};
            
            # we don't accept opts from user, since we handle them all ourselves
            # and they must not be overridden
            my $plot_opts = '-r ' . $self->inputs->{fasta_gc_stats_file}->[0]->path;
            
            foreach my $s_file (@{ $self->inputs->{stats_files} }) {
                # we need to know some related info on this stats file;
                # find it in the graph database under the vrtrack schema or die
                my $vrtrack = VRPipe::Schema->create('VRTrack');
                my $vrstats_file = $vrtrack->get_file($s_file->protocolless_path, $s_file->protocol);
                $self->throw($s_file->path . " was not in the graph database") unless $vrstats_file;
                my $vr_lane = $vrstats_file->closest('VRTrack', 'Lane', direction => 'incoming');
                my ($vr_stats) = $vrstats_file->related(outgoing => { type => 'summary_stats' });
                my $prefix = $vr_lane->unique() || $s_file->basename;
                
                my @vr_plot_params;
                $self->output_file(temporary => 1, basename => $prefix . '-quals.gp', type => 'txt');
                my $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-quals.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Qualities' });
                $self->output_file(temporary => 1, basename => $prefix . '-quals2.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-quals2.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Qualities' });
                $self->output_file(temporary => 1, basename => $prefix . '-quals3.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-quals3.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Qualities' });
                
                $self->output_file(temporary => 1, basename => $prefix . '-quals-hm.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-quals-hm.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Qualities' });
                
                if ($vr_stats->properties->{'insert size average'} > 0) {
                    $self->output_file(temporary => 1, basename => $prefix . '-insert-size.gp', type => 'txt');
                    $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-insert-size.png', type => 'bin');
                    push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Insert Size' });
                }
                
                $self->output_file(temporary => 1, basename => $prefix . '-gc-content.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-gc-content.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'GC Content' });
                
                $self->output_file(temporary => 1, basename => $prefix . '-gc-depth.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-gc-depth.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'GC Depth' });
                
                $self->output_file(temporary => 1, basename => $prefix . '-acgt-cycles.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-acgt-cycles.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'ACGT Cycles' });
                
                $self->output_file(temporary => 1, basename => $prefix . '-coverage.gp', type => 'txt');
                $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-coverage.png', type => 'bin');
                push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Coverage' });
                
                if ($vr_stats->properties->{'number of insertions'} > 0 && $vr_stats->properties->{'number of deletions'} > 0) {
                    $self->output_file(temporary => 1, basename => $prefix . '-indel-dist.gp', type => 'txt');
                    $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-indel-dist.png', type => 'bin');
                    push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Indel distribution' });
                    $self->output_file(temporary => 1, basename => $prefix . '-indel-cycles.gp', type => 'txt');
                    $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-indel-cycles.png', type => 'bin');
                    push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Indels per cycle' });
                }
                
                if ($vr_stats->properties->{'options'} =~ /(?<!\S)-r(?!\S)/) {
                    $self->output_file(temporary => 1, basename => $prefix . '-mism-per-cycle.gp', type => 'txt');
                    $ofile = $self->output_file(output_key => 'bamstats_plots', basename => $prefix . '-mism-per-cycle.png', type => 'bin');
                    push(@vr_plot_params, { path => $ofile->path->stringify, type => 'png', caption => 'Mismatches per cycle' });
                }
                
                $self->output_file(temporary => 1, basename => $prefix . '.html', type => 'txt');
                
                # first remove any old plots still attached to the stats file
                my @existing_plots = $vrstats_file->related(outgoing => { type => 'bamstats_plot' });
                foreach my $plot (@existing_plots) {
                    $vrstats_file->divorce_from($plot);
                }
                
                foreach my $params (@vr_plot_params) {
                    my $path = delete $params->{path};
                    $self->relate_input_to_output($vrstats_file, 'bamstats_plot', $path, $params);
                }
                
                my $cat = $s_file->cat_cmd;
                my $cmd = qq[$cat | $plot_bamstats $plot_opts -p $prefix -];
                $self->dispatch([$cmd, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bamstats_plots => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'png files produced by plot-bamstats, with a caption in the metadata',
                min_files   => 8,
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Uses plot-bamstats to draw plots from the associated samtools stats file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
