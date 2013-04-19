
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
    
    method submit_args (VRPipe::Requirements :$requirements!, Str|File :$stdo_file!, Str|File :$stde_file!, Str :$cmd!, VRPipe::PersistentArray :$array?) {
        # access the requirements object and build up the string based on
        # memory, time, cpu etc.
        my $instance_type = $self->determine_queue($requirements);
        # *** ...
        my $megabytes          = $requirements->memory;
        my $requirments_string = "instance $instance_type memory $megabytes";
        my $cpus               = $requirements->cpus;
        if ($cpus > 1) {
            $requirments_string .= " cpus $cpus";
        }
        
        return qq[$requirments_string cmd '$cmd'];
    }
    
    sub submit {
        my ($self, %args) = @_;
        my $instance_type = $args{instance} || $self->throw("No instance supplied");
        my $megabytes     = $args{memory}   || $self->throw("No memory supplied");
        my $cmd           = $args{cmd}      || $self->throw("No cmd supplied");
        
        #*** not yet implemented
        warn "would submit cmd [$cmd] to instance [$instance_type], requiring [$megabytes]MB\n";
        
        # is there already an instance running with 'room' for our command?
        my $instance;
        my @current_instances = $ec2->describe_instances({
                'image-id'          => $ami,
                'availability-zone' => $availability_zone,
                'instance-type'     => $instance_type,
                #'tag:Key'                        => 'Value'
            }
        );
        
        foreach my $possible (@current_instances) {
            #*** decide if this instance can be used
            
            last if $instance;
        }
        
        unless ($instance) {
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
        }
        
        if ($instance) {
            my $instance_id = $instance->instanceId;
            warn "selected instance $instance_id\n";
            
            #*** not yet implemented ...
            
            # print "Job x is submitted\n";
            
            sleep(20);
            warn "will terminate the instance\n";
            $instance->terminate;
        }
        else {
            $self->throw("Could not find or create an instance to submit the command to");
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
    
    method switch_queue (PositiveInt $sid, Str $new_queue) {
        # we don't support queue switching
        $self->throw("Queue Switching is not supported (and should not be necessary) for the ec2 scheduler");
    }
    
    method get_scheduler_id {
        #*** not yet implemented
        return -1;
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        # we don't have any concept of a job 'array', so don't deal with indexes
        return 1;
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
    
    method kill_sids (ArrayRef $sid_aids) {
        #*** not yet implemented
        
        my @sids;
        foreach my $sid_aid (@$sid_aids) {
            my ($sid, $aid) = @$sid_aid;
            my $id = $aid ? qq{"$sid\[$aid\]"} : $sid;
            push(@sids, $id);
            
            if (@sids == 500) {
                #system("bkill @sids");
                @sids = ();
            }
        }
        
        #system("bkill @sids") if @sids;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        #*** not yet implemented
        
        my $id = $sid; # we don't support Job arrays, so ignore $aid
        
        my $status;
        return $status || 'UNKNOWN';
    }
    
    method run_time (PositiveInt $sid, Int $aid) {
        my $id = $sid;
        
        my ($start_epoch, $end_epoch);
        
        #*** not yet implemented
        
        $start_epoch || return 0;
        $end_epoch ||= time();
        return $end_epoch - $start_epoch;
    }
    
    method command_status (Str :$cmd, PositiveInt :$max?) {
        my $count            = 0;
        my @running_sid_aids = ();
        my @to_kill;
        my $job_name_prefix = $self->_job_name($cmd);
        
        #*** not yet implemented
        
        if (@to_kill) {
            $self->kill_sids(@to_kill);
        }
        
        return ($count, \@running_sid_aids);
    }
}

1;
