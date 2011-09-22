use VRPipe::Base;

class VRPipe::Schedulers::lsf with VRPipe::SchedulerMethodsRole {
    method start_command {
        return 'bjobs'; #*** actually, I've no idea how to start lsf
    }
    method stop_command {
        return 'bjobs'; #*** actually, I've no idea how to stop lsf
    }
    
    method submit_command {
        return 'bsub';
    }
    
    method submit_args (VRPipe::Requirements :$requirements, File :$stdo_file, File :$stde_file, Str :$cmd, VRPipe::PersistentArray :$array?) {
        # access the requirments object and build up the string based on memory,
        # time, cpu etc.
        my $queue = $self->determine_queue($requirements);
        # *** ...
        my $megabytes = $requirements->memory;
        my $m = $megabytes * 1000;
        my $requirments_string = "-q $queue -M$m -R 'select[mem>$megabytes] rusage[mem=$megabytes]'";
        
        # work out the scheduler output locations and how to pass on the
        # scheduler array index to the perl cmd
        my $index_spec = '';
        my $array_def = '';
        my $output_string;
        if ($array) {
            $index_spec = ''; #*** something that gives the index to be shifted into perl -e; for LSF we leave it empty and will pick up an env var elsewhere instead
            $output_string = "-o $stdo_file.\%I -e $stde_file.\%I";
            my $size = $array->size;
            my $uniquer = $array->id;
            $array_def = qq{-J "vrpipe$uniquer\[1-$size]" };
        }
        else {
            $output_string = "-o $stdo_file -e $stde_file";
        }
        
        #*** could solve issue of having to have _aid and _hid in Submission
        #    purely to allow us to find where the lsf output went by having
        #    output go to named pipe that stores directly in db against the
        #    submission id. Same could be done for Job output. Would it work
        #    across the farm?
        
        return qq[$array_def$output_string $requirments_string '$cmd$index_spec'];
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        #*** sanger specific, need a way of having this module ask its queue
        #    related questions during siteconfig setup, so we can avoid hard
        #    coding here
        return $requirements->time > 12 ? 'long' : 'normal';
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        $index ? return $index : return $ENV{LSB_JOBINDEX};
    }
    
    method get_sid (Str $cmd) {
        my $output = `$cmd`;
        my ($sid) = $output =~ /Job \<(\d+)\> is submitted/;
        if ($sid) {
            return $sid;
        }
        else {
            $self->throw("Failed to submit to scheduler");
        }
    }
    
    method kill_sid (PositiveInt $sid, Int $aid, PositiveInt $secs = 30) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid;
        my $t = time();
        while (1) {
            last if time() - $t > $secs;
            
            #*** fork and kill child if over time limit?
            my $status = $self->sid_status($sid, $aid);
            last if ($status eq 'UNKNOWN' || $status eq 'DONE' || $status eq 'EXIT');
            
            system("bkill $id");
            
            sleep(1);
        }
        return 1;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid; # when aid is 0, it was not a job array
        open(my $bfh, "bjobs $id |") || $self->warn("Could not call bjobs $id");
        my $status;
        if ($bfh) {
            while (<$bfh>) {
                if (/^$sid\s+\S+\s+(\S+)/) {
                    $status = $1;
                }
            }
            close($bfh);
        }
        
        return $status || 'UNKNOWN'; # *** needs to return a word in a defined vocabulary suitable for all schedulers
    }
}

1;