use VRPipe::Base;

class VRPipe::Schedulers::local with VRPipe::SchedulerMethodsRole {
    method submit_command {
        return 'vrpipe-local_scheduler submit';
    }
    
    method submit_args (VRPipe::Requirements :$requirements, File :$stdo_file, File :$stde_file, Str :$cmd, VRPipe::PersistentArray :$array?) {
        my $array_def = '';
        my $output_string;
        if ($array) {
            $output_string = "-o $stdo_file.\%I -e $stde_file.\%I";
            my $size = $array->size;
            $array_def = "-a $size";
        }
        else {
            $output_string = "-o $stdo_file -e $stde_file";
        }
        
        return qq[$array_def$output_string '$cmd'];
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        return 'local';
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        $index ? return $index : return $ENV{VRPIPE_LOCAL_JOBINDEX};
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
            
            system("vrpipe-local_scheduler kill $id");
            
            sleep(1);
        }
        return 1;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        my $id = $aid ? qq{"$sid\[$aid\]"} : $sid; # when aid is 0, it was not a job array
        open(my $bfh, "vrpipe-local_scheduler jobs $id |") || $self->warn("Could not call vrpipe-local_scheduler jobs $id");
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