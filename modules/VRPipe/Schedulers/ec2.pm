
=head1 NAME

VRPipe::Schedulers::ec2 - interface to ec2

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This class provides L<Amazon
EC2|http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/Welcome.html>-specific
command lines for use by L<VRPipe::Scheduler>.

It depends upon VM::EC2, which you must manually install.

For this to work the server must be running on an AMI that boots up to an
environment with a working VRPipe installation and all the software you need to
run. The VRPipe installation must be configured to use ec2 as the job
scheduler, and you must have provided the access and secret keys (which can be
found at https://portal.aws.amazon.com/gp/aws/securityCredentials?) during the
'perl Build.PL' phase of VRPipe installation.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Schedulers::ec2 with VRPipe::SchedulerMethodsRole {
    # eval this so that test suite can pass syntax check on this module when
    # VM::EC2 is not installed
    eval "use VM::EC2;";
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::Interface::BackEnd;
    use Net::SSH qw(ssh ssh_cmd);
    use POSIX qw(ceil);
    
    #*** are instance type details not query-able? Do we have to hard-code it?
    our %instance_types = (
        'm1.small'   => [1, 1700,  1,    0.065, 0.065], # cores, MB, ECU(speed), cost/hr, cost/hr/speed
        'm1.medium'  => [1, 3750,  2,    0.130, 0.07],
        'm1.large'   => [2, 7500,  2,    0.260, 0.13],
        'm1.xlarge'  => [4, 15000, 2,    0.520, 0.26],
        'm3.xlarge'  => [4, 15000, 3.25, 0.550, 0.17],
        'm3.2xlarge' => [8, 30000, 3.25, 1.100, 0.34],
        't1.micro'   => [1, 613,   1,    0.020, 0.02],
        'm2.xlarge'  => [2, 17100, 3.25, 0.460, 0.14],
        'm2.2xlarge' => [4, 34200, 3.25, 0.920, 0.28],
        'm2.4xlarge' => [8, 68400, 3.25, 1.840, 0.57],
        'c1.medium'  => [2, 1700,  2.5,  0.165, 0.07],
        'c1.xlarge'  => [8, 7000,  2.5,  0.660, 0.26]
    );
    
    # we expect that the majority of what we run will be single cpu, cpu
    # intensive jobs. Therefore we want to pick the type that will get the work
    # done quickest considering only one of its cores, at the lowest cost:
    # cost/hr/speed. But when the cost/hr/speed is very close for types that are
    # very different in speed, it makes more sense to pick the faster one since
    # fewer hours may be used. Because of this we just hard-code a preferred
    # order that makes the most sense.
    our @ordered_types = ('t1.micro', 'c1.medium', 'm1.small', 'm1.medium', 'm1.large', 'm2.xlarge', 'm3.xlarge', 'c1.xlarge', 'm1.xlarge', 'm2.2xlarge', 'm3.2xlarge', 'm2.4xlarge');
    
    our %queues;
    our $access_key        = $vrp_config->ec2_access_key;
    our $secret_key        = $vrp_config->ec2_secret_key;
    our $url               = $vrp_config->ec2_url;
    our ($region)          = $url =~ /ec2\.(.+?)\.amazonaws/;
    our $key_name          = $vrp_config->ec2_private_key_name;
    our $ec2               = VM::EC2->new(-access_key => $access_key, -secret_key => $secret_key, -region => $region);
    our $meta              = $ec2->instance_metadata;
    our $ami               = $meta->imageId;
    our @security_groups   = $meta->securityGroups;
    our $availability_zone = $meta->availabilityZone;
    our $deployment        = VRPipe::Persistent::SchemaBase->database_deployment;
    our $backend           = VRPipe::Interface::BackEnd->new(deployment => $deployment);
    
    method start_command {
        return 'sleep 1'; #*** not really applicable
    }
    
    method stop_command {
        return 'sleep 1'; #*** not really applicable
    }
    
    method submit_command {
        # we call a method in this module to submit
        my $perl = VRPipe::Interface::CmdLine->vrpipe_perl_command($deployment);
        return $perl . ' -MVRPipe::Schedulers::ec2 -e "VRPipe::Schedulers::ec2->submit(@ARGV)"';
    }
    
    method submit_args (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, PositiveInt :$count = 1) {
        # access the requirements object and build up the string based on
        # memory and cpu (other reqs do not apply)
        my $instance_type      = $self->determine_queue($requirements);
        my $megabytes          = $requirements->memory;
        my $requirments_string = "instance $instance_type memory $megabytes";
        my $cpus               = $requirements->cpus;
        if ($cpus > 1) {
            $requirments_string .= " cpus $cpus";
        }
        
        return qq[$requirments_string count $count cmd '$cmd'];
    }
    
    sub submit {
        my ($self, %args) = @_;
        my $instance_type = $args{instance} || $self->throw("No instance supplied");
        my $megabytes     = $args{memory}   || $self->throw("No memory supplied");
        my $count         = $args{count}    || $self->throw("No count supplied");
        my $cmd           = $args{cmd}      || $self->throw("No cmd supplied");
        my $cpus          = $args{cpus}     || 1;
        
        warn "will submit cmd [$cmd] to instance [$instance_type] $count times, requiring [$megabytes]MB\n";
        
        #*** there are surely optimisations and caching that could be done here...
        for (1 .. $count) {
            # is there already an instance running with 'room' for our command?
            my $instance;
            my @current_instances = $ec2->describe_instances({
                    'image-id'            => $ami,
                    'availability-zone'   => $availability_zone,
                    'instance-type'       => $instance_type,
                    'instance-state-name' => 'running'
                }
            );
            
            my $available_cpus   = $instance_types{$instance_type}->[0];
            my $available_memory = $instance_types{$instance_type}->[1];
            my $own_pdn          = $meta->privateDnsName;
            foreach my $possible (@current_instances) {
                # see what vrpipe-handler processes are running on this instance
                # (searching our own job table for jobs running on this host isn't
                #  good enough, since the handler may not have started running a job
                #  yet)
                my $pdn = $possible->privateDnsName;
                next if $pdn eq $own_pdn; # submit() will be called by the server, and we don't want any handlers running on the same instance as the server
                my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
                warn "will search for processes running on $host\n";
                my $processes = ssh_with_return($host, qq[ps xj | grep vrpipe-handler]) || '';
                my %pgids;
                foreach my $process (split("\n", $processes)) {
                    next if $process =~ /grep/;
                    my ($pgid) = $process =~ /\s*\d+\s+\d+\s+(\d+)/;
                    my ($r)    = $process =~ /-r (\d+) /;
                    $pgids{$pgid} = $r || 0;
                }
                
                my $cpus_used   = 0;
                my $memory_used = 0;
                while (my ($pgid, $rid) = each %pgids) {
                    my $req = VRPipe::Requirements->get(id => $rid) if $rid;
                    $cpus_used += $req ? $req->cpus : $available_cpus;
                    last if $cpus_used >= $available_cpus;
                    
                    # get the total memory used by all processes in this process
                    # group
                    my $processes = ssh_with_return($host, qq[ps xj | grep $pgid]) || '';
                    my $this_memory_used = 0;
                    foreach my $process (split("\n", $processes)) {
                        next if $process =~ /grep/;
                        my ($pid, $this_pgid) = $process =~ /\s*\d+\s+(\d+)\s+(\d+)/;
                        next unless $this_pgid == $pgid;
                        my $grep = ssh_with_return($host, qq[grep VmRSS /proc/$pid/status 2>/dev/null]);
                        my $grep_bytes;
                        if ($grep && $grep =~ /(\d+) kB/) {
                            $grep_bytes = $1 * 1024;
                            $this_memory_used += ceil($grep_bytes / 1048576);
                        }
                    }
                    
                    my $this_memory_reserved = $req ? $req->memory : $this_memory_used;
                    if ($this_memory_used > $this_memory_reserved) {
                        warn "will try to kill pgid $pgid because it is using too much memory\n";
                        ssh($pdn, qq[kill -TERM -$pgid]);
                        $memory_used += $this_memory_used;
                    }
                    else {
                        $memory_used += $this_memory_reserved;
                    }
                    last if $memory_used >= $available_memory;
                }
                
                warn "$host had used $cpus_used/$available_cpus cpus and $memory_used/$available_memory memory\n";
                if ($cpus <= $available_cpus - $cpus_used && $megabytes <= $available_memory - $memory_used) {
                    $instance = $possible;
                    warn "will use $host\n";
                    last;
                }
            }
            
            unless ($instance) {
                warn "no suitable instances, will spawn a new one\n";
                # launch a new instance
                ($instance) = $ec2->run_instances(
                    -image_id               => $ami,
                    -instance_type          => $instance_type,
                    -client_token           => $ec2->token,
                    -key_name               => $key_name,
                    -security_group         => \@security_groups,
                    -availability_zone      => $availability_zone,
                    -min_count              => 1,
                    -max_count              => 1,
                    -termination_protection => 0,
                    -shutdown_behavior      => 'terminate'
                ) or $self->throw($ec2->error_str);
                
                $ec2->wait_for_instances($instance);
                my $status = $instance->current_status;
                $self->throw("Created a new instance but it didn't start running normally") unless $status eq 'running';
                warn "started up instance ", $instance->instanceId, " which has host ", $instance->privateDnsName, "\n";
            }
            
            if ($instance) {
                my $instance_id = $instance->instanceId;
                my $ip          = $instance->privateIpAddress;
                warn "selected instance $instance_id at $ip\n";
                
                #*** can't figure out how to both detatch and return immediately
                # from running the command over ssh, and get back the pid of the
                # cmd we just started, so we'll have to note the pids already on
                # the machine so we can detect afterwards what new pid was created
                my $processes = ssh_with_return($ip, qq[ps xj | grep vrpipe-handler]) || '';
                my %existing_pgids;
                foreach my $process (split("\n", $processes)) {
                    next if $process =~ /grep/;
                    next unless $process =~ /$cmd/;
                    my ($pgid) = $process =~ /\s*\d+\s+\d+\s+(\d+)/;
                    $existing_pgids{$pgid} = 1;
                }
                
                ssh($ip, qq[sh -c "( ( nohup $cmd &>/dev/null ) & )"]);
                
                $processes = ssh_with_return($ip, qq[ps xj | grep vrpipe-handler]) || '';
                my $pgid;
                foreach my $process (split("\n", $processes)) {
                    next if $process =~ /grep/;
                    next unless $process =~ /$cmd/;
                    my ($this_pgid) = $process =~ /\s*\d+\s+\d+\s+(\d+)/;
                    unless (exists $existing_pgids{$this_pgid}) {
                        $pgid = $this_pgid;
                        last;
                    }
                }
                
                if ($pgid) {
                    print "Job <$ip:$pgid> is submitted\n";
                }
                else {
                    $self->throw("Failed to launch cmd on $ip via ssh");
                }
            }
            else {
                $self->throw("Could not find or create an instance to submit the command to");
            }
        }
    }
    
    method terminate_old_instances (PositiveInt $max_do_nothing_time = 3600) {
        warn "will check for instances that can be terminated\n";
        my @all_instances = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'running'
            }
        );
        my $own_pdn = $meta->privateDnsName;
        foreach my $instance (@all_instances) {
            my $pdn = $instance->privateDnsName;
            next if $pdn eq $own_pdn; # don't terminate ourselves - the server that calls this method won't have any handlers running on it
            my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
            my $jobs = VRPipe::Job->search({ host => $host, heartbeat => { '>=' => DateTime->from_epoch(epoch => time() - $max_do_nothing_time) } });
            next if $jobs;
            warn "will terminate instance $host\n";
            $instance->terminate;
        }
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        # based on the requirements we want to select an appropriate ec2
        # instance type to run the job on
        my $megabytes = $requirements->memory;
        my $cpus      = $requirements->cpus;
        
        # we hard-coded a preferred order for types above; select the first
        # type in the list that meets our requirements
        foreach my $type (@ordered_types) {
            my ($available_cpus, $available_megabytes) = @{ $instance_types{$type} };
            next if $available_cpus < $cpus;
            next if $available_megabytes < $megabytes;
            return $type;
        }
        
        $self->throw("No EC2 instance type is compatible with running jobs requiring $cpus cpus and $megabytes MB of memory");
    }
    
    method queue_time (VRPipe::Requirements $requirements) {
        # we can run for an unlimited time on all instance types
        return 31536000;
    }
    
    method switch_queue (Str $sid, Str $new_queue) {
        # we don't support queue switching
        $self->throw("Queue Switching is not supported (and should not be necessary) for the ec2 scheduler");
    }
    
    method get_scheduler_id {
        my $ip   = $meta->privateIpAddress;
        my $pgid = getpgrp(0);
        return "$ip:$pgid";
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        # we don't have any concept of a job 'array', so don't deal with indexes
        return 0;
    }
    
    method get_sid (Str $cmd) {
        my $output = `$cmd`;
        my ($sid) = $output =~ /Job \<(.+?)\> is submitted/;
        
        if ($sid) {
            return $sid;
        }
        else {
            $self->throw("Failed to submit to scheduler");
        }
    }
    
    method kill_sids (ArrayRef $sid_aids) {
        foreach my $sid_aid (@$sid_aids) {
            my ($sid) = @$sid_aid;
            my ($ip, $pgid) = split(':', $sid);
            ssh($ip, qq[kill -TERM -$pgid]);
        }
    }
    
    method sid_status (Str $sid, Int $aid) {
        my ($ip, $pgid) = split(':', $sid);
        
        my $processes = ssh_with_return($ip, qq[ps xj | grep vrpipe-handler]) || '';
        my $found = 0;
        foreach my $process (split("\n", $processes)) {
            my ($this_pgid) = $process =~ /\s*\d+\s+\d+\s+(\d+)/;
            if ($this_pgid == $pgid) {
                $found = 1;
                last;
            }
        }
        
        return $found ? 'RUN' : 'UNKNOWN';
    }
    
    method run_time (Str $sid, Int $aid) {
        my ($ip, $pgid) = split(':', $sid);
        
        my $processes = ssh_with_return($ip, qq[ps xj | grep vrpipe-handler]) || '';
        my $pid = 0;
        foreach my $process (split("\n", $processes)) {
            my ($this_pid, $this_pgid) = $process =~ /\s*\d+\s+(\d+)\s+(\d+)/;
            if ($this_pgid == $pgid) {
                $pid = $this_pid;
                last;
            }
        }
        $pid || return 0;
        
        my $etime = ssh_with_return($ip, qq[ps -p $pid -o etime=]);
        # [[dd-]hh:]mm:ss
        my ($d, $h, $m, $s) = $etime =~ /(?:(?:(\d+)-)?(\d+):)?(\d+):(\d+)/;
        $d ||= 0;
        $h ||= 0;
        my $seconds = ($d * 24 * 60 * 60) + ($h * 60 * 60) + ($m * 60) + $s;
        
        return $seconds;
    }
    
    method command_status (Str :$cmd, PositiveInt :$max?) {
        # we have no concept of a pending job, since our submit() method just
        # immediately starts running cmds on a node: we don't have to care about
        # $max or killing things here
        
        #*** there must be some optimisation involving caching or tags or redis
        # that can be done here?...
        
        # look through all running instances and see what's running on them
        my $count            = 0;
        my @running_sid_aids = ();
        my @all_instances    = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'running'
            }
        );
        foreach my $instance (@all_instances) {
            my $ip = $instance->privateIpAddress;
            
            my $processes = ssh_with_return($ip, qq[ps xj | grep vrpipe-handler]) || '';
            my %pgids;
            foreach my $process (split("\n", $processes)) {
                next unless $process =~ /$cmd/;
                my ($pgid) = $process =~ /\s*\d+\s+\d+\s+(\d+)/;
                $pgids{$pgid} = 1;
            }
            
            $count += keys %pgids;
            foreach my $pgid (keys %pgids) {
                push(@running_sid_aids, ["$ip:$pgid", 0]);
            }
        }
        
        return ($count, \@running_sid_aids);
    }
    
    # we wrap Net::SSH's ssh_cmd method to avoid complications with our tied
    # STDERR
    sub ssh_with_return {
        my ($host, $cmd) = @_;
        my $tied = tied *STDERR ? 1 : 0;
        untie *STDERR if $tied;
        my $return = ssh_cmd($host, $cmd);
        $backend->log_stderr() if $tied;
        return $return;
    }
}

1;
