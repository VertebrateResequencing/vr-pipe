
=head1 NAME

VRPipe::Schedulers::sge_ec2 - interface to SGE on EC2

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This class provides L<Sub (Son of)
GridEngine|https://arc.liv.ac.uk/trac/SGE>-specific command lines for use by
L<VRPipe::Scheduler>.

It differs from the sge scheduler by launching or terminating EC2 instances and
adding/removing them from the SGE cluster according to current demand.

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

class VRPipe::Schedulers::sge_ec2 extends VRPipe::Schedulers::sge {
    our $ec2_scheduler = VRPipe::SchedulerMethodsFactory->create('ec2', {});
    our $ec2 = $VRPipe::Schedulers::ec2::ec2;
    
    sub periodic_method {
        return 'launch_extra_instances_and_terminate_old_instances';
    }
    
    sub on_exit_method {
        return 'terminate_all_instances';
    }
    
    method launch_extra_instances_and_terminate_old_instances (Str :$deployment!) {
        ## launch extra instances
        
        # look at what we have queuing to get a count of how many of each req
        # we need to run, then launch the appropriate instances for those reqs
        # and install SGE on them
        my %launched_hosts;
        
        ## terminate old instances
        
        # start by seeing which hosts are not currently running any jobs
        open(my $qfh, 'qhost -j -ncb |') || $self->throw("Unable to open pipe from qhost -j -ncb");
        <$qfh>;
        <$qhf>;
        <$qfh>; # ignore the first 3 lines
        my %hosts;
        my $host;
        while (<$qfh>) {
            if (/^\S+/) {
                my @cols = split;
                $host = $cols[0];
                $hosts{$host}->{load} = $cols[3];
            }
            elsif (/^\s+\d+\s+/) {
                $hosts{$host}->{jobs}++;
            }
        }
        close($qfh) || $self->throw("Unable to close pipe from qhost -j -ncb");
        
        # consider empty hosts to be those that also haven't run a job in the
        # past hour
        my $max_do_nothing_time = $deployment eq 'production' ? 3600 : 300;
        my @empty_hosts;
        while (my ($host, $details) = each %hosts) {
            # ones we just now launched probably won't be running anything yet
            next if exists $launched_hosts{$host};
            
            # 'load' is a value that drops off over time once nothing is
            # running on the machine; in production wait until it as the lowest
            # value of 0.01
            if ($deployment eq 'production') {
                next if $details->{load} > 0.01;
            }
            
            # obviously it's not empty if it is running a job right now
            next if $details->{jobs};
            
            # see if a job has run on this host in the past hour
            next if VRPipe::Job->search({ host => $host, heartbeat => { '>=' => DateTime->from_epoch(epoch => time() - $max_do_nothing_time) } });
            
            # this host has been sitting there doing nothing for ages
            push(@empty_hosts, $host);
        }
        
        # stop SGE on the empty hosts
        my $sge_config = $sge_config = $self->_alter_config($sge_config, 'EXEC_HOST_LIST_RM', join(' ', @empty_hosts));
        #[edit /gluster/vrpipe/.sge_config to set host in EXEC_HOST_LIST_RM]
        system('cd $SGE_ROOT; ./inst_sge -ux -auto ' . $sge_config) && $self->throw("Failed to run: inst_sge -ux -auto $sge_config");
        unlink($sge_config);
        
        # terminate the instances
        my $own_pdn = $meta->privateDnsName;
        foreach my $host (@empty_hosts) {
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
}

1;
