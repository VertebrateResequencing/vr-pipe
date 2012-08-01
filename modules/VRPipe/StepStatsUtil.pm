
=head1 NAME

VRPipe::StepStatsUtil - analyse past performance of Steps

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Typically a L<VRPipe::Step> will define some command line(s) that needs to be
executed. For efficiency and good job scheduling reasons, B<VRPipe> would like
to know how much time and memory these command lines usually require to
execute. StepStatsUtil looks at the past time and memory used (as recorded in
L<VRPipe::StepStat> objects) and works out appropriate values to use in the
future.

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

class VRPipe::StepStatsUtil {
    use POSIX qw(ceil);
    
    our %means;
    our %percentiles;
    
    has 'step' => (is       => 'rw',
                   isa      => 'VRPipe::Step',
                   required => 1);
    
    method _percentile (Str $column! where { $_ eq 'memory' || $_ eq 'time' }, Int $percent! where { $_ > 0 && $_ < 100 }, VRPipe::PipelineSetup $pipelinesetup?) {
        # did we already work this out in the past hour?
        my $step      = $self->step;
        my $time      = time();
        my $store_key = $step->id . $column . $percent . ($pipelinesetup ? $pipelinesetup->id : 0);
        my $p_data    = $percentiles{$store_key};
        if ($p_data && $p_data->[0] + 3600 > $time) {
            return ($p_data->[1], $p_data->[2]);
        }
        
        my @search_args = (step => $step->id, $pipelinesetup ? (pipelinesetup => $pipelinesetup->id) : ());
        my ($count, $percentile) = (0, 0);
        $count = VRPipe::StepStats->search({@search_args});
        if ($count) {
            ($percentile) = VRPipe::StepStats->get_column_values($column, {@search_args}, { order_by => { -desc => [$column] }, rows => 1, offset => sprintf("%0.0f", ($count / 100) * (100 - $percent)) });
        }
        
        $percentiles{$store_key} = [$time, $count, $percentile] if $count > 500;
        return ($count, $percentile);
    }
    
    method percentile_seconds (Int :$percent!, VRPipe::PipelineSetup :$pipelinesetup?) {
        return $self->_percentile('time', $percent, $pipelinesetup ? ($pipelinesetup) : ());
    }
    
    method percentile_memory (Int :$percent!, VRPipe::PipelineSetup :$pipelinesetup?) {
        return $self->_percentile('memory', $percent, $pipelinesetup ? ($pipelinesetup) : ());
    }
    
    method _mean (Str $column! where { $_ eq 'memory' || $_ eq 'time' }, VRPipe::PipelineSetup $pipelinesetup?) {
        # did we already work this out in the past hour?
        my $step      = $self->step;
        my $time      = time();
        my $store_key = $step->id . $column . ($pipelinesetup ? $pipelinesetup->id : 0);
        my $mean_data = $means{$store_key};
        if ($mean_data && $mean_data->[0] + 3600 > $time) {
            return ($mean_data->[1], $mean_data->[2], $mean_data->[3]);
        }
        
        # get the mean and sd using little memory
        my ($count, $mean, $sd) = (0, 0, 0);
        my $pager = VRPipe::StepStats->get_column_values_paged($column, { step => $step->id, $pipelinesetup ? (pipelinesetup => $pipelinesetup->id) : () });
        while (my $stats = $pager->next) {
            foreach my $stat (@$stats) {
                $count++;
                if ($count == 1) {
                    $mean = $stat;
                    $sd   = 0;
                }
                else {
                    my $old_mean = $mean;
                    $mean += ($stat - $old_mean) / $count;
                    $sd += ($stat - $old_mean) * ($stat - $mean);
                }
            }
        }
        
        if ($count) {
            $mean = sprintf("%0.0f", $mean);
            $sd   = sprintf("%0.0f", sqrt($sd / $count));
        }
        
        $means{$store_key} = [$time, $count, $mean, $sd] if $count > 500;
        return ($count, $mean, $sd);
    }
    
    method mean_seconds (VRPipe::PipelineSetup :$pipelinesetup?) {
        return $self->_mean('time', $pipelinesetup ? ($pipelinesetup) : ());
    }
    
    method mean_memory (VRPipe::PipelineSetup :$pipelinesetup?) {
        return $self->_mean('memory', $pipelinesetup ? ($pipelinesetup) : ());
    }
    
    method _recommended (Str $method! where { $_ eq 'memory' || $_ eq 'time' }, VRPipe::PipelineSetup $pipelinesetup?) {
        # if we've seen enough previous results, recommend 95th percentile
        # rounded up to nearest 100
        my ($count, $percentile) = $self->_percentile($method, 95, $pipelinesetup ? ($pipelinesetup) : ());
        if ($count >= 3) {
            if ($percentile % 100) {
                return ceil($percentile / 100) * 100;
            }
            else {
                return $percentile;
            }
        }
        return;
    }
    
    method recommended_memory (VRPipe::PipelineSetup :$pipelinesetup?) {
        my $mem = $self->_recommended('memory', $pipelinesetup ? ($pipelinesetup) : ()) || return;
        
        # recommend at least 100MB though
        if ($mem < 100) {
            $mem = 100;
        }
        return $mem;
    }
    
    method recommended_time (VRPipe::PipelineSetup :$pipelinesetup?) {
        my $seconds = $self->_recommended('time', $pipelinesetup ? ($pipelinesetup) : ()) || return;
        
        # convert to hrs, rounded up
        return ceil($seconds / 60 / 60);
    }
}

1;
