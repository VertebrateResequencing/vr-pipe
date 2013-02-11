
=head1 NAME

VRPipe::Schedulers::lsf - interface to LSF

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This class provides L<Amazon
EC2|http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/Welcome.html>-specific
command lines for use by L<VRPipe::Scheduler>.

It depends upon VM::EC2, which you must manually install.

For this to work you must have exactly one AMI with the name 'vrpipe-install',
and it must be private to you. You can create such an AMI by starting with our
public AMI and customising for your needs. The AMI must boot up to an
environment with a working VRPipe installation and all the software you need to
run. The VRPipe installation must be configured to use ec2 as the job
scheduler, and you must have provided the access and secret keys (which can be
found at https://portal.aws.amazon.com/gp/aws/securityCredentials?).

You must also have a security group called 'vrpipe-security' which at a minimum
must allow ssh between ec2 instances.

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
    
    our %queues;
    our $ami;
    our $access_key = $vrp_config->ec2_access_key;
    our $secret_key = $vrp_config->ec2_secret_key;
    our $url        = $vrp_config->ec2_url;
    our ($region)   = $url =~ /ec2\.(.+?)\.amazonaws/;
    our $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
    
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
        my $instance_type = $args{'instance'} || $self->throw("No instance supplied");
        my $megabytes     = $args{'memory'}   || $self->throw("No memory supplied");
        my $cmd           = $args{'cmd'}      || $self->throw("No cmd supplied");
        
        #*** not yet implemented
        warn "would submit cmd [$cmd] to instance [$instance_type], requiring [$megabytes]MB\n";
    }
    
    method determine_queue (VRPipe::Requirements $requirements) {
        # based on the requirements we want to select an appropriate ec2
        # instance type to run the job on
        
        #my $megabytes = $requirements->memory;
        #my $chosen_queue;
        #*** not yet implemented
        #return $chosen_queue;
        
        #*** for now we just hard-code m1.medium, which is 1 core, 3.7GB
        
        return 'm1.medium';
    }
    
    method queue_time (VRPipe::Requirements $requirements) {
        # we can run for an unlimited time on all instance types
        return 31536000;
    }
    
    method switch_queue (PositiveInt $sid, Str $new_queue) {
        # we don't support queue switching
        $self->throw("Queue Switching is not supported (and should not be necessary) for the ec2 scheduler");
    }
    
    method get_1based_index (Maybe[PositiveInt] $index?) {
        # we don't have any concept of a job 'array', so don't deal with indexes
        return 1;
    }
    
    method get_sid (Str $cmd) {
        my $ec2;
        unless ($ami) {
            # get user's private vrpipe-install AMI if we don't already know it
            $ec2 = VM::EC2->new(-access_key => $access_key, -secret_key => $secret_key, -region => $region);
            my @amis = $ec2->describe_images(-owner => "self", -filter => { name => "vrpipe-install" }) or $self->throw("desribe_images call failed: " . $ec2->error_str);
            unless (@amis == 1) {
                $self->throw("There was not exactly 1 vrpipe-install AMI owned by you");
            }
            $ami = $amis[0]->imageId;
            warn "got ami $ami\n";
        }
        
        #*** not yet implemented
        my $sid;
        
        #my $output = `$cmd`;
        #my ($sid) = $output =~ /Job \<(\d+)\> is submitted/;
        
        if ($sid) {
            return $sid;
        }
        else {
            $self->throw("Failed to submit to scheduler");
        }
    }
    
    method kill_sid (PositiveInt $sid, Int $aid, PositiveInt $secs = 30) {
        #*** not yet implemented
        
        my $t = time();
        while (1) {
            last if time() - $t > $secs;
            
            #*** fork and kill child if over time limit?
            my $status = $self->sid_status($sid, $aid);
            last if ($status eq 'UNKNOWN' || $status eq 'DONE' || $status eq 'EXIT');
            
            #system("bkill $id");
            
            sleep(1);
        }
        return 1;
    }
    
    method all_status {
        my %status = ();
        
        #*** not yet implemented
        # $status{$sid} = 'RUN';
        
        return %status;
    }
    
    method sid_status (PositiveInt $sid, Int $aid) {
        my $status;
        #*** not yet implemented
        return $status || 'UNKNOWN'; # *** needs to return a word in a defined vocabulary suitable for all schedulers
    }
}

1;
