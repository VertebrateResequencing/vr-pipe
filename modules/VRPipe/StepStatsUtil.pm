use VRPipe::Base;

class VRPipe::StepStatsUtil {
    use POSIX;
    
    our %means;
    
    has 'step' => (is => 'rw',
                   isa => 'VRPipe::Step',
                   required => 1);
    
    method _mean (Str $column! where { $_ eq 'memory' || $_ eq 'time' }, VRPipe::PipelineSetup $pipelinesetup?) {
        # did we already work this out in the past hour?
        my $step = $self->step;
        my $time = time();
        my $store_key = $step->id.$column. ($pipelinesetup ? $pipelinesetup->id : 0);
        my $mean_data = $means{$store_key};
        if ($mean_data && $mean_data->[0] + 3600 > $time) {
            return ($mean_data->[1], $mean_data->[2], $mean_data->[3]);
        }
        
        # get the mean and sd using little memory
        my ($count, $mean, $sd) = (0, 0, 0);
        my $schema = $step->result_source->schema;
        my $rs = $schema->resultset('StepStats')->search({ step => $step->id, $pipelinesetup ? (pipelinesetup => $pipelinesetup->id) : () });
        my $rs_column = $rs->get_column($column);
        while (my $stat = $rs_column->next) { # using $rs_column instead of $rs is >60x faster with 100k+ rows
            $count++;
            if ($count == 1) {
                $mean = $stat;
                $sd = 0;
            }
            else {
                my $old_mean = $mean;
                $mean += ($stat - $old_mean) / $count;
                $sd += ($stat - $old_mean) * ($stat - $mean);
            }
        }
        
        if ($count) {
            $mean = sprintf("%0.0f", $mean);
            $sd = sprintf("%0.0f", sqrt($sd / $count));
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
        my ($count, $mean, $sd) = $self->_mean($method, $pipelinesetup ? ($pipelinesetup) : ());
        
        # if we've seen enough completed submissions, recommend mean + 2sd
        if ($count >= 3) {
            return $mean + 2 * $sd;
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