
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

class VRPipe::Schedulers::ec2 extends VRPipe::Schedulers::local {
    # eval this so that test suite can pass syntax check on this module when
    # VM::EC2 is not installed
    eval "use VM::EC2;";
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use DateTime::Format::Natural;
    
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
    
    sub submit {
        my ($self, %args) = @_;
        my $queue     = $args{queue}  || $self->throw("No queue supplied");
        my $megabytes = $args{memory} || $self->throw("No memory supplied");
        my $count     = $args{count}  || $self->throw("No count supplied");
        my $cmd       = $args{cmd}    || $self->throw("No cmd supplied");
        my $cpus      = $args{cpus}   || 1;
        my $cwd       = $args{cwd};
        
        # we do not maintain state or a queue or anything like that. We just
        # check if it is possible to run the $cmd, and do so up to $count times.
        # If we fail to run the command $count times, we don't care; vrpipe-
        # server will just resubmit in the future
        
        #*** there are surely optimisations and caching that could be done
        # here...
        
        # do a pass to see if there's room to run the cmd on any of our ips
        my $usable_ips = $self->_usable_ips($queue, $megabytes, $cpus, $count) || return;
        
        # spawn new instances if needed
        unless (@$usable_ips == $count) {
            my $needed = $count - @$usable_ips;
            warn "insufficient suitable instances, will spawn $needed new ones\n";
            # by default people are limited to a max of 20 instances:
            # http://www.phacai.com/increase-ec2-instance-quota
            # so we have to handle a possible error here
            my %run_instance_args = (
                -image_id               => $ami,
                -instance_type          => $queue,
                -client_token           => $ec2->token,
                -key_name               => $key_name,
                -security_group         => \@security_groups,
                -availability_zone      => $availability_zone,
                -min_count              => $needed,
                -max_count              => $needed,
                -termination_protection => 0,
                -shutdown_behavior      => 'terminate'
            );
            my @new_instances = $ec2->run_instances(%run_instance_args);
            
            unless (@new_instances) {
                my $error = $ec2->error_str;
                if ($error =~ /instances exceeds your current quota of (\d+)/) {
                    my $max = $1;
                    my @all_instances = $ec2->describe_instances({ 'instance-state-name' => 'running' });
                    $count = $max - @all_instances;
                    if ($count == 0) {
                        $backend->log("Unable to spawn any new instances; consider increasing your quota: http://aws.amazon.com/contact-us/ec2-request/");
                    }
                    else {
                        $backend->log("Your EC2 account has an instance quota of $max; we need $needed more instances, but will launch $count new ones instead");
                        $run_instance_args{'-min_count'} = $count;
                        $run_instance_args{'-max_count'} = $count;
                        @new_instances                   = $ec2->run_instances(%run_instance_args);
                        unless (@new_instances) {
                            $backend->log("Failed to launch $count new instances: " . $ec2->error_str);
                        }
                    }
                }
                else {
                    $backend->log("Failed to launch $needed new instances: " . $error);
                }
            }
            
            $ec2->wait_for_instances(@new_instances) if @new_instances;
            
            # check they're all fine
            my $available_cpus   = $instance_types{$queue}->[0];
            my $available_memory = $instance_types{$queue}->[1];
            foreach my $instance (@new_instances) {
                my $iip    = $instance->privateIpAddress;
                my $status = $instance->current_status;
                unless ($status eq 'running') {
                    $backend->log("Created a new ec2 instance at $iip but it didn't start running normally");
                    next;
                }
                warn "started up instance at $iip\n";
                
                # wait for it to become responsive to ssh
                my $max_tries  = 240;
                my $responsive = 0;
                for (1 .. $max_tries) {
                    my $return = $backend->ssh($iip, 'echo ssh_working');
                    if ($return && $return =~ /ssh_working/) {
                        $responsive = 1;
                        warn "the instance was responsive to ssh\n";
                        last;
                    }
                    sleep(1);
                }
                unless ($responsive) {
                    $backend->log("Newly launched instance $iip is not responding to ssh");
                    next;
                }
                
                my $cpus_used   = 0;
                my $memory_used = 0;
                while (@$usable_ips < $count) {
                    if ($cpus <= $available_cpus - $cpus_used && $megabytes <= $available_memory - $memory_used) {
                        push(@$usable_ips, $iip);
                        $cpus_used   += $cpus;
                        $memory_used += $megabytes;
                    }
                    else {
                        last;
                    }
                }
            }
        }
        
        # start a handler on each usable ip (_usable_ips would have returned the
        # same ip more than once if that up could support running more than 1
        # handler)
        $self->_start_new_handlers($usable_ips, $cmd, $cwd ? ($cwd) : ());
    }
    
    method _cluster_ips (Str $queue?) {
        my @current_instances = $ec2->describe_instances({
                'image-id'          => $ami,
                'availability-zone' => $availability_zone,
                $queue ? ('instance-type' => $queue) : (),
                'instance-state-name' => 'running'
            }
        );
        
        my @ips;
        my $own_ip = $meta->privateIpAddress;
        foreach my $instance (@current_instances) {
            my $this_ip = $instance->privateIpAddress;
            next if $own_ip eq $this_ip; # submit() will be called by the server, and we don't want any handlers running on the same instance as the server
            push(@ips, $this_ip);
        }
        return @ips;
    }
    
    method terminate_old_instances (Str :$deployment!) {
        my $dt_parser = DateTime::Format::Natural->new;
        
        warn "will check for instances that can be terminated\n";
        my $max_do_nothing_time = $deployment eq 'production' ? 3600 : 300;
        my @all_instances = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'running'
            }
        );
        my $own_pdn = $meta->privateDnsName;
        foreach my $instance (@all_instances) {
            # don't terminate ourselves - the server that calls this method
            # won't have any handlers running on it
            my $pdn = $instance->privateDnsName;
            next if $pdn eq $own_pdn;
            
            my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
            
            # don't terminate if we only just now spawned it and maybe are still
            # waiting for it to become responsive to ssh
            my $dstr = $instance->launchTime;
            $dstr =~ s/(\d)T(\d)/$1 $2/;
            $dstr =~ s/\.\d+Z$//;
            my $dt      = $dt_parser->parse_datetime($dstr);
            my $elapsed = time() - $dt->epoch;
            next if $elapsed < 300;
            
            # don't terminate an instance that has a handler running on it right
            # now
            next if $self->_handler_processes($host);
            
            # don't terminate if the instance has recently run a Job
            next if VRPipe::Job->search({ host => $host, heartbeat => { '>=' => DateTime->from_epoch(epoch => time() - $max_do_nothing_time) } });
            
            warn "will terminate instance $host\n";
            $instance->terminate;
        }
    }
    
    method terminate_all_instances {
        warn "will find all instances to terminate them\n";
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
            $instance->terminate;
        }
    }
    
    sub periodic_method {
        return 'terminate_old_instances';
    }
    
    sub on_exit_method {
        return 'terminate_all_instances';
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
    
    method get_scheduler_id {
        my $ip   = $meta->privateIpAddress;
        my $pgid = getpgrp(0);
        return "$ip:$pgid";
    }
}

1;
