
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
    use DateTime;
    use DateTime::Format::Natural;
    use LWP::UserAgent;
    use JSON::XS;
    use VRPipe::Persistent::InMemory;
    
    #*** are instance type details not query-able? Do we have to hard-code it?
    # the costs seem stable, but differs by region; I only include them here to
    # help me rank them manually
    our %instance_types = (
        't1.micro'    => [1,  613,    1,    0.020, 0.02], # cores, MB, ECU(speed), cost/hr, cost/hr/speed
        'm1.small'    => [1,  1700,   1,    0.065, 0.065],
        'm1.medium'   => [1,  3750,   2,    0.130, 0.07],
        'm1.large'    => [2,  7500,   2,    0.260, 0.13],
        'm1.xlarge'   => [4,  15000,  2,    0.520, 0.26],
        'm2.xlarge'   => [2,  17100,  3.25, 0.460, 0.14],
        'm2.2xlarge'  => [4,  34200,  3.25, 0.920, 0.28],
        'm2.4xlarge'  => [8,  68400,  3.25, 1.840, 0.57],
        'm3.xlarge'   => [4,  15000,  3.25, 0.550, 0.17],
        'm3.2xlarge'  => [8,  30000,  3.25, 1.100, 0.34],
        'c1.medium'   => [2,  1700,   2.5,  0.165, 0.07],
        'c1.xlarge'   => [8,  7000,   2.5,  0.660, 0.26],
        'cc1.4xlarge' => [8,  23000,  4.19, 1.300, 0.31],
        'cc2.8xlarge' => [16, 60500,  5.5,  2.700, 0.49],
        'cr1.8xlarge' => [16, 244000, 5.5,  3.750, 0.68],
        'cg1.4xlarge' => [8,  22500,  4.19, 2.360, 0.56],
        'hi1.4xlarge' => [8,  60500,  4.37, 3.410, 0.78],
        'hs1.8xlarge' => [16, 117000, 2.19, 4.900, 2.24]
    );
    
    our %price_type_to_instance_type = (
        uODI   => { u => 't1.micro' },
        stdODI => {
            sm  => 'm1.small',
            med => 'm1.medium',
            lg  => 'm1.large',
            xl  => 'm1.xlarge'
        },
        hiMemODI => {
            xl    => 'm2.xlarge',
            xxl   => 'm2.2xlarge',
            xxxxl => 'm2.4xlarge'
        },
        secgenstdODI => {
            xl  => 'm3.xlarge',
            xxl => 'm3.2xlarge'
        },
        hiCPUODI => {
            med => 'c1.medium',
            xl  => 'c1.xlarge'
        },
        clusterHiMemODI => { xxxxxxxxl => 'cr1.8xlarge' },
        clusterComputeI => {
            xxxxl     => 'cc1.4xlarge',
            xxxxxxxxl => 'cc2.8xlarge'
        },
        clusterGPUI => { xxxxl     => 'cg1.4xlarge' },
        hiIoODI     => { xxxxl     => 'hi1.4xlarge' },
        hiStoreODI  => { xxxxxxxxl => 'hs1.8xlarge' }
    );
    
    our %price_region_to_area = (
        'us-east'    => 'us-east-1',
        'us-west-2'  => 'us-west-2',
        'us-west'    => 'us-west-1',
        'eu-ireland' => 'eu-west-1',
        'apac-sin'   => 'ap-southeast-1',
        'apac-syd'   => 'ap-southeast-2',
        'apac-tokyo' => 'ap-northeast-1',
        'sa-east-1'  => 'sa-east-1'
    );
    
    our $access_key         = $vrp_config->ec2_access_key;
    our $secret_key         = $vrp_config->ec2_secret_key;
    our $url                = $vrp_config->ec2_url;
    our ($region)           = $url =~ /ec2\.(.+?)\.amazonaws/;
    our $key_name           = $vrp_config->ec2_private_key_name;
    our $max_instances      = $vrp_config->ec2_max_instances;
    our $spot_price_percent = $vrp_config->ec2_spot_price_percent;
    our $spot_force_percent = $vrp_config->ec2_spot_force_percent;
    our @ordered_types      = split(',', $vrp_config->ec2_instance_types);
    our $ec2                = VM::EC2->new(-access_key => $access_key, -secret_key => $secret_key, -region => $region);
    our $meta               = $ec2->instance_metadata;
    our $ami                = $meta->imageId;
    our @security_groups    = $meta->securityGroups;
    our $availability_zone  = $meta->availabilityZone;
    our $deployment         = VRPipe::Persistent::SchemaBase->database_deployment;
    our $backend            = VRPipe::Interface::BackEnd->new(deployment => $deployment);
    
    our $initialized = 0;
    our %instance_on_demand_price;
    our $price_list_url = 'http://aws.amazon.com/ec2/pricing/pricing-on-demand-instances.json';
    our %instance_lowest_spot_price;
    our $lowest_spot_price_check_time;
    
    my $im = VRPipe::Persistent::InMemory->new();
    
    around initialize_for_server {
        return if $initialized;
        
        if ($spot_price_percent) {
            $im->log("[ec2scheduler] Spot requests enabled");
            
            # grab the latest on-demand instance pricing for our region
            my $area = $region;
            $area =~ s/[a-z]$//;
            
            my $ua       = LWP::UserAgent->new;
            my $response = $ua->get($price_list_url);
            if ($response->is_success) {
                my $prices_json = $response->decoded_content;
                my $json        = JSON::XS->new->utf8->canonical;
                my $hash        = $json->decode($prices_json);
                my $region_data = $hash->{config}->{regions};
                
                foreach my $rd (@$region_data) {
                    my $this_area = $price_region_to_area{ $rd->{region} };
                    $this_area && $this_area eq $area || next;
                    
                    foreach my $type_ref (@{ $rd->{instanceTypes} }) {
                        my $base_type = $type_ref->{type};
                        
                        foreach my $size_ref (@{ $type_ref->{sizes} }) {
                            my $size = $size_ref->{size};
                            my $type = $price_type_to_instance_type{$base_type}->{$size} || next;
                            
                            foreach my $vc (@{ $size_ref->{valueColumns} }) {
                                $vc->{name} eq "linux" || next;
                                my $usd = $vc->{prices}->{USD};
                                next unless $usd && $usd =~ /^[\d\.]+$/; # can be N/A
                                $instance_on_demand_price{$type} = $usd;
                                $im->debug("[ec2scheduler] on-demand price for $type = $usd");
                            }
                        }
                    }
                    
                    last;
                }
            }
            else {
                $im->log("[ec2scheduler] Spot requests disabled because I could not download the on-demand price list at $price_list_url :" . $response->status_line);
                $spot_price_percent = 0;
            }
            
            $self->_get_lowest_spot_prices;
        }
        else {
            $im->log("[ec2scheduler] Spot requests disabled");
        }
        
        # call parent's initialize method
        $self->$orig();
        
        $initialized = 1;
    }
    
    method _get_lowest_spot_prices {
        return unless $spot_price_percent;
        return if $spot_force_percent;
        
        # only do this every 24hrs
        my $t = time();
        return unless (!$lowest_spot_price_check_time || $t > ($lowest_spot_price_check_time + 86400));
        
        # find the lowest spot prices in the last 24hrs, so we know if the
        # user's $spot_price_percent of on-demand price should be used, or if we
        # should use the lowest spot (which ever is higher)
        $im->debug("[ec2scheduler] Will check the lowest spot prices over the last 24hrs");
        
        my $yesterday = DateTime->from_epoch(epoch => time() - 86400);
        while (my ($type, $odp) = each %instance_on_demand_price) {
            my @spot_price_history = $ec2->describe_spot_price_history(-start_time => "$yesterday", -instance_type => $type, -availability_zone => $availability_zone, -filter => { "product-description" => "Linux/UNIX*" });
            
            my $min_price = $odp * 2;
            foreach my $h (@spot_price_history) {
                my $spot_price = $h->spot_price;
                if ($spot_price < $min_price) {
                    $min_price = $spot_price;
                }
            }
            
            $im->debug("[ec2scheduler] lowest spot price for $type = $min_price");
            $instance_lowest_spot_price{$type} = $min_price;
        }
        
        $lowest_spot_price_check_time = $t;
    }
    
    sub submit {
        my ($self, %args) = @_;
        my $type      = $args{queue}  || $self->throw("No queue supplied");
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
        my $usable_ips = $self->_usable_ips($type, $megabytes, $cpus, $count) || return;
        
        # spawn new instances if needed
        unless (@$usable_ips == $count) {
            my $needed = $count - @$usable_ips;
            if ($needed > 0) {
                my @new_instances = $self->launch_instances($type, $needed);
                $count = $count - ($needed - @new_instances);
                my $available_cpus   = $instance_types{$type}->[0];
                my $available_memory = $instance_types{$type}->[1];
                foreach my $instance (@new_instances) {
                    my $iip         = $instance->privateIpAddress;
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
        }
        
        # start a handler on each usable ip (_usable_ips would have returned the
        # same ip more than once if that up could support running more than 1
        # handler)
        $self->_start_new_handlers($usable_ips, $cmd, $cwd ? ($cwd) : ());
    }
    
    method launch_instances (Str $type, PositiveInt $needed) {
        # we won't launch more instances than the user-set max
        my @pending_instances = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'pending'
            }
        );
        my @running_instances = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'running'
            }
        );
        my $current_instances = @pending_instances + @running_instances;
        my $allowed_instances = $max_instances - $current_instances + 1; # +1 because one of them will be the instance we're launching from
        
        $im->debug("[ec2scheduler] Got a request to launch $needed $type instances; allowed to launch $allowed_instances");
        if ($needed > $allowed_instances) {
            $needed = $allowed_instances;
        }
        return if $needed <= 0;
        
        my %run_instance_args = (
            -image_id               => $ami,
            -instance_type          => $type,
            -key_name               => $key_name,
            -security_group         => \@security_groups,
            -availability_zone      => $availability_zone,
            -min_count              => $needed,
            -max_count              => $needed,
            -termination_protection => 0,
            -shutdown_behavior      => 'terminate'
        );
        
        my @new_instances;
        if ($spot_price_percent) {
            my $odp = $instance_on_demand_price{$type};
            $odp || $self->throw("[ec2scheduler] No on-demand price found for instance type $type");
            my $bid = $spot_price_percent == 100 ? $odp : sprintf("%.3f", (($odp / 100) * $spot_price_percent));
            unless ($spot_force_percent) {
                my $min = $instance_lowest_spot_price{$type};
                if ($min > $bid) {
                    $im->debug("[ec2scheduler] Minimum historical spot price was greater than $spot_price_percent\% of on-demand price ($bid), so will bid at $min");
                    $bid = $min;
                }
            }
            
            if ($bid > $odp) {
                # if the bid is greater than on-demand price, switch to
                # on-demand mode
                $im->debug("[ec2scheduler] Bid $bid is greater than on-demand price; switching to on-demand mode");
                @new_instances = $ec2->run_instances(%run_instance_args);
            }
            else {
                # make a spot request and wait until it is fulfilled
                # (ideally we would not wait since in the time it takes to
                # fulfil we may no longer need this many instances, but we don't
                # have enough info in this call to know our needs for all
                # instance types so can't properly manage previous requests)
                my %spot_args = (%run_instance_args, -spot_price => $bid, -instance_count => $needed);
                my @requests = $ec2->request_spot_instances(%spot_args);
                
                if ($ec2->is_error) {
                    # if we're over the instance count limit, keep retrying
                    # until it works
                    my $error = $ec2->error;
                    if ($error->code eq "MaxSpotInstanceCountExceeded") {
                        my @pending_instances = $ec2->describe_instances({
                                'availability-zone'   => $availability_zone,
                                'instance-state-name' => 'pending'
                            }
                        );
                        my @running_instances = $ec2->describe_instances({
                                'availability-zone'   => $availability_zone,
                                'instance-state-name' => 'running'
                            }
                        );
                        my $total_instances  = @pending_instances + @running_instances;
                        my @pending_requests = $ec2->describe_spot_instance_requests({ state => "open", 'launched-availability-zone' => $availability_zone });
                        my @active_requests  = $ec2->describe_spot_instance_requests({ state => "active", 'launched-availability-zone' => $availability_zone });
                        my $total_requests   = @pending_requests + @active_requests;
                        
                        $spot_args{'-instance_count'}--;
                        while ($spot_args{'-instance_count'} > 0) {
                            @requests = $ec2->request_spot_instances(%spot_args);
                            
                            if ($ec2->is_error) {
                                if ($error->code eq "MaxSpotInstanceCountExceeded") {
                                    $spot_args{'-instance_count'}--;
                                }
                                else {
                                    $im->log("[ec2scheduler] Failed to request $needed new spot instances: " . $ec2->error_str);
                                    last;
                                }
                            }
                            else {
                                # try and work out what the quota must be and
                                # then set $max_instances so that we can avoid
                                # this code path in future
                                my $quota = $spot_args{'-instance_count'} + $total_requests;
                                $max_instances = $quota + ($total_instances - $total_requests);
                                
                                $im->log("[ec2scheduler] Failed to launch $needed more instances due to exceeding your quota of ~$quota, but launched " . $spot_args{'-instance_count'} . " new ones instead; in future I will only try to launch up to $max_instances instances");
                                
                                last;
                            }
                        }
                    }
                    else {
                        $im->log("[ec2scheduler] Failed to request $needed new spot instances: " . $ec2->error_str);
                    }
                }
                
                foreach my $request (@requests) {
                    while (1) {
                        my $status = $request->current_status;
                        my $state  = $request->state;
                        
                        if ($state ne 'open' && $state ne 'active') {
                            $im->log("[ec2scheduler] Spot instance request id: " . $request->spotInstanceRequestId . " entered unexpected state '$state'; skipping this request");
                            last;
                        }
                        
                        $im->debug("[ec2scheduler] spot instance request id: " . $request->spotInstanceRequestId . " | state: $state; status: $status");
                        
                        my $instanceid = $request->instanceId;
                        if ($instanceid) {
                            $im->debug("[ec2scheduler] Spot instance request " . $request->spotInstanceRequestId . " resulted in instance $instanceid");
                            push(@new_instances, $request->instance);
                            last;
                        }
                        sleep(5);
                    }
                }
            }
        }
        else {
            @new_instances = $ec2->run_instances(%run_instance_args);
        }
        
        $im->debug("[ec2scheduler] Tried to launch $needed $type instances; actually launched " . scalar(@new_instances));
        
        # by default people are limited to a max of 20 instances:
        # http://www.phacai.com/increase-ec2-instance-quota
        # so we have to handle a possible error here. Unfortunately I don't
        # know of a way of finding out the quota, and the error messages change
        # on some random basis, with a different meaning for the numbers in the
        # error message. So if it's a quota issue we'll just lower our requested
        # number until it works
        if ($ec2->is_error) {
            if ($ec2->error->code eq 'InstanceLimitExceeded') {
                my $total_instances;
                unless ($spot_price_percent) {
                    my @pending_instances = $ec2->describe_instances({
                            'availability-zone'   => $availability_zone,
                            'instance-state-name' => 'pending'
                        }
                    );
                    my @running_instances = $ec2->describe_instances({
                            'availability-zone'   => $availability_zone,
                            'instance-state-name' => 'running'
                        }
                    );
                    $total_instances = @pending_instances + @running_instances;
                }
                
                my $count = $needed - 1;
                while ($count > 0) {
                    $run_instance_args{'-min_count'} = $count;
                    $run_instance_args{'-max_count'} = $count;
                    @new_instances                   = $ec2->run_instances(%run_instance_args);
                    if ($ec2->is_error) {
                        if ($ec2->error->code eq 'InstanceLimitExceeded') {
                            $count--;
                            $im->debug("[ec2scheduler] got InstanceLimitExceeded, will try $count instances");
                        }
                        else {
                            $im->log("[ec2scheduler] Failed to launch $needed more instances due to exceeding your quota, and also failed to launch $count new instances: " . $ec2->error_str);
                            last;
                        }
                    }
                    else {
                        unless ($spot_price_percent) {
                            # try and work out what the quota must be and then
                            # set $max_instances so that we can avoid this code
                            # path in future
                            my $quota = $count + $total_instances;
                            $max_instances = $quota - ($total_instances - $current_instances + 1);
                            
                            $im->log("[ec2scheduler] Failed to launch $needed more instances due to exceeding your quota of ~$quota, but launched $count new ones instead; in future I will only try to launch up to $max_instances instances");
                        }
                        
                        last;
                    }
                }
                if ($count == 0) {
                    $im->log("[ec2scheduler] Failed to launch $needed more instances due to exceeding your quota; consider increasing your quota: http://aws.amazon.com/contact-us/ec2-request/");
                }
            }
            else {
                $im->log("[ec2scheduler] Failed to launch $needed new instances: " . $ec2->error_str);
            }
        }
        
        # check they're all fine
        my @good_instances;
        if (@new_instances) {
            $ec2->wait_for_instances(@new_instances);
            
            foreach my $instance (@new_instances) {
                my $iip    = $instance->privateIpAddress;
                my $status = $instance->current_status;
                unless ($status eq 'running') {
                    $im->log("[ec2scheduler] Created a new ec2 instance at $iip but it didn't start running normally: will terminate it");
                    $instance->terminate;
                    next;
                }
                
                # lock the instance so we don't try and terminate it before it
                # becomes responsive to SSH
                my $lock_key = 'starting_instance.' . $iip;
                $im->note($lock_key, forget_after => 360);
                
                # wait for it to become responsive to ssh
                my $max_time   = time() + 300;
                my $responsive = 0;
                while (time() < $max_time) {
                    my $return;
                    eval {
                        local $SIG{ALRM} = sub { die "alarm\n" };
                        alarm 15;
                        $return = $backend->ssh($iip, 'echo ssh_working');
                        alarm 0;
                    };
                    
                    if ($return && $return =~ /ssh_working/) {
                        $responsive = 1;
                        last;
                    }
                    sleep(1);
                }
                unless ($responsive) {
                    $im->log("[ec2scheduler] Newly launched instance $iip is not responding to ssh: will terminate it");
                    $instance->terminate;
                    $im->forget_note($lock_key);
                    next;
                }
                
                $im->log("[ec2scheduler] Started up instance at $iip");
                push(@good_instances, $instance);
                
                # update the lock timeout so that we don't terminate an instance
                # within a minute of starting it
                $im->note($lock_key, forget_after => 60);
            }
        }
        
        return @good_instances;
    }
    
    method _cluster_ips (Str $type?) {
        my @current_instances;
        unless ($type) {
            @current_instances = $ec2->describe_instances({
                    'image-id'            => $ami,
                    'availability-zone'   => $availability_zone,
                    'instance-state-name' => 'running'
                }
            );
        }
        else {
            # $type is an instance type, but don't just get back instances
            # that exactly match that type; also get instances that match or
            # exceed the specs of that type
            my @itypes;
            my ($needed_cores, $needed_mem) = @{ $instance_types{$type} };
            foreach my $itype (@ordered_types) {
                next if $itype eq $type;
                my ($available_cores, $available_mem) = @{ $instance_types{$itype} };
                next if $available_cores < $needed_cores;
                next if $available_mem < $needed_mem;
                push(@itypes, $itype);
            }
            
            foreach my $itype ($type, @itypes) {
                push(
                    @current_instances,
                    $ec2->describe_instances({
                            'image-id'            => $ami,
                            'availability-zone'   => $availability_zone,
                            'instance-type'       => $itype,
                            'instance-state-name' => 'running'
                        }
                    )
                );
            }
        }
        
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
        
        my $max_do_nothing_time = $deployment eq 'production' ? 3600 : 300;
        my @all_instances = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'running'
            }
        );
        my $own_pdn = $meta->privateDnsName;
        my ($own_host) = $own_pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
        foreach my $instance (@all_instances) {
            # don't terminate ourselves - the server that calls this method
            # won't have any handlers running on it
            my $pdn = $instance->privateDnsName;
            my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
            next if $host eq $own_host;
            
            # don't terminate if we only just now spawned it and maybe are still
            # waiting for it to become responsive to ssh
            next if $im->noted('starting_instance.' . $instance->privateIpAddress);
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
            next if VRPipe::Job->search({ host => $host, end_time => { '>=' => DateTime->from_epoch(epoch => time() - $max_do_nothing_time) } });
            
            $self->log("[ec2scheduler] Will terminate instance $host");
            $instance->terminate;
        }
        
        # this is our periodic method, so we'll ask the vrpipe-server process
        # to run _get_lowest_spot_prices
        return ('_get_lowest_spot_prices');
    }
    
    method terminate_all_instances {
        my @all_instances = $ec2->describe_instances({
                'image-id'            => $ami,
                'availability-zone'   => $availability_zone,
                'instance-state-name' => 'running'
            }
        );
        my $own_pdn = $meta->privateDnsName;
        my ($own_host) = $own_pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
        foreach my $instance (@all_instances) {
            # don't terminate ourselves - the server that calls this method
            # won't have any handlers running on it
            my $pdn = $instance->privateDnsName;
            my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
            next if $host eq $own_host;
            
            # don't terminate an instance that has a handler running on it right
            # now - possibly spawned by a production server
            next if $self->_handler_processes($host);
            
            $self->log("[ec2scheduler] Will terminate instance $host");
            $instance->terminate;
        }
    }
    
    sub periodic_method {
        return 'terminate_old_instances';
    }
    
    sub on_exit_method {
        return 'terminate_all_instances';
    }
    
    method determine_queue (VRPipe::Requirements $requirements, Int $global_max = 0) {
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
