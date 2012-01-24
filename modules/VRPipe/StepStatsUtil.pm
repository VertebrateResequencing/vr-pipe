use VRPipe::Base;

class VRPipe::StepStatsUtil {
    use POSIX;
    
    has 'step' => (is => 'rw',
                   isa => 'VRPipe::Step',
                   required => 1);
    
    method stepstats (VRPipe::PipelineSetup :$pipelinesetup?) {
        my $rs = $self->_stepstats_rs($pipelinesetup ? ($pipelinesetup) : ());
        my @step_stats;
        while (my $ss = $rs->next) {
            push(@step_stats, $ss);
        }
        return @step_stats;
    }
    method _stepstats_rs (VRPipe::PipelineSetup $pipelinesetup?) {
        my $step = $self->step;
        my $schema = $step->result_source->schema;
        my $rs = $schema->resultset('StepStats')->search({ step => $step->id, $pipelinesetup ? (pipelinesetup => $pipelinesetup->id) : () });
        return $rs;
    }
    
    method _mean (Str $method! where { $_ eq 'memory' || $_ eq 'time' }, VRPipe::PipelineSetup $pipelinesetup?) {
        # get the mean and sd using little memory
        my ($count, $mean, $sd) = (0, 0, 0);
        my $rs = $self->_stepstats_rs($pipelinesetup ? ($pipelinesetup) : ());
        while (my $ss = $rs->next) {
            my $stat = $ss->$method() || next;
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
        
        # if we've seen enough completed submissions, recommend mean + 2sd, but
        # at least 100MB
        if ($count >= 3) {
            my $rec = $mean + 2 * $sd;
            if ($rec < 100) {
                $rec = 100;
            }
            return $rec;
        }
        return;
    }
    method recommended_memory (VRPipe::PipelineSetup :$pipelinesetup?) {
        return $self->_recommended('memory', $pipelinesetup ? ($pipelinesetup) : ());
    }
    method recommended_time (VRPipe::PipelineSetup :$pipelinesetup?) {
        my $seconds = $self->_recommended('time', $pipelinesetup ? ($pipelinesetup) : ()) || return;
        
        # convert to hrs, rounded up
        return ceil($seconds / 60 / 60);
    }
}

1;