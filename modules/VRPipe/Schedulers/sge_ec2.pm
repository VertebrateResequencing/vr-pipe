
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

Note that it alters configuration of your SGE installation, changing things
like your complex attributes, parallel environment, queues and groups. (And so
the user that runs vrpipe-server must have permission to do these things.)

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
    use VRPipe::Config;
    use Path::Class;
    use DateTime;
    use POSIX qw(ceil);
    use VRPipe::Persistent::InMemory;
    
    our $ec2_scheduler = VRPipe::SchedulerMethodsFactory->create('ec2', {});
    no warnings 'once';
    our $ec2 = $VRPipe::Schedulers::ec2::ec2;
    my $im = VRPipe::Persistent::InMemory->new();
    
    my $vrp_config      = VRPipe::Config->new();
    my $sge_config_file = $vrp_config->sge_config_file;
    my $log_dir_method  = VRPipe::Persistent::SchemaBase->database_deployment . '_logging_directory';
    our $sge_confs_dir = dir($vrp_config->$log_dir_method(), '.sge_confs');
    our $complex_attribs = <<COMPLEX;
#name               shortcut   type      relop requestable consumable default  urgency
#--------------------------------------------------------------------------------------
arch                a          STRING    ==    YES         NO         NONE     0
calendar            c          STRING    ==    YES         NO         NONE     0
cpu                 cpu        DOUBLE    >=    YES         NO         0        0
display_win_gui     dwg        BOOL      ==    YES         NO         0        0
h_core              h_core     MEMORY    <=    YES         NO         0        0
h_cpu               h_cpu      TIME      <=    YES         NO         0:0:0    0
h_data              h_data     MEMORY    <=    YES         NO         0        0
h_fsize             h_fsize    MEMORY    <=    YES         NO         0        0
h_rss               h_rss      MEMORY    <=    YES         NO         0        0
h_rt                h_rt       TIME      <=    YES         NO         0:0:0    0
h_stack             h_stack    MEMORY    <=    YES         NO         0        0
h_vmem              h_vmem     MEMORY    <=    YES         NO         0        0
hostname            h          HOST      ==    YES         NO         NONE     0
load_avg            la         DOUBLE    >=    NO          NO         0        0
load_long           ll         DOUBLE    >=    NO          NO         0        0
load_medium         lm         DOUBLE    >=    NO          NO         0        0
load_short          ls         DOUBLE    >=    NO          NO         0        0
m_core              core       INT       <=    YES         NO         0        0
m_socket            socket     INT       <=    YES         NO         0        0
m_thread            thread     INT       <=    YES         NO         0        0
m_topology          topo       STRING    ==    YES         NO         NONE     0
m_topology_inuse    utopo      STRING    ==    YES         NO         NONE     0
mem_free            mf         MEMORY    <=    FORCED      YES        0        0
mem_total           mt         MEMORY    <=    YES         NO         0        0
mem_used            mu         MEMORY    >=    YES         NO         0        0
min_cpu_interval    mci        TIME      <=    NO          NO         0:0:0    0
np_load_avg         nla        DOUBLE    >=    NO          NO         0        0
np_load_long        nll        DOUBLE    >=    NO          NO         0        0
np_load_medium      nlm        DOUBLE    >=    NO          NO         0        0
np_load_short       nls        DOUBLE    >=    NO          NO         0        0
num_proc            p          INT       ==    YES         NO         0        0
qname               q          STRING    ==    YES         NO         NONE     0
rerun               re         BOOL      ==    NO          NO         0        0
s_core              s_core     MEMORY    <=    YES         NO         0        0
s_cpu               s_cpu      TIME      <=    YES         NO         0:0:0    0
s_data              s_data     MEMORY    <=    YES         NO         0        0
s_fsize             s_fsize    MEMORY    <=    YES         NO         0        0
s_rss               s_rss      MEMORY    <=    YES         NO         0        0
s_rt                s_rt       TIME      <=    YES         NO         0:0:0    0
s_stack             s_stack    MEMORY    <=    YES         NO         0        0
s_vmem              s_vmem     MEMORY    <=    YES         NO         0        0
seq_no              seq        INT       ==    NO          NO         0        0
slots               s          INT       <=    YES         YES        1        1000
swap_free           sf         MEMORY    <=    YES         NO         0        0
swap_rate           sr         MEMORY    >=    YES         NO         0        0
swap_rsvd           srsv       MEMORY    >=    YES         NO         0        0
swap_total          st         MEMORY    <=    YES         NO         0        0
swap_used           su         MEMORY    >=    YES         NO         0        0
tmpdir              tmp        STRING    ==    NO          NO         NONE     0
virtual_free        vf         MEMORY    <=    YES         YES        0        0
virtual_total       vt         MEMORY    <=    YES         NO         0        0
virtual_used        vu         MEMORY    >=    YES         NO         0        0
COMPLEX
    our $queue_defaults = <<QUEUE;
seq_no                0
load_thresholds       np_load_avg=1.75
suspend_thresholds    NONE
nsuspend              1
suspend_interval      00:05:00
priority              0
min_cpu_interval      00:05:00
processors            UNDEFINED
qtype                 BATCH INTERACTIVE
ckpt_list             NONE
pe_list               threaded
rerun                 FALSE
tmpdir                /tmp
shell                 /bin/sh
prolog                NONE
epilog                NONE
shell_start_mode      posix_compliant
starter_method        NONE
suspend_method        NONE
resume_method         NONE
terminate_method      NONE
notify                00:00:60
owner_list            NONE
user_lists            NONE
xuser_lists           NONE
subordinate_list      NONE
projects              NONE
xprojects             NONE
calendar              NONE
initial_state         default
s_rt                  INFINITY
h_rt                  INFINITY
s_cpu                 INFINITY
h_cpu                 INFINITY
s_fsize               INFINITY
h_fsize               INFINITY
s_stack               INFINITY
h_stack               INFINITY
s_core                INFINITY
h_core                INFINITY
QUEUE
    our $scheduling = <<SCHED;
algorithm                         default
schedule_interval                 0:0:15
maxujobs                          0
queue_sort_method                 load
job_load_adjustments              NONE
load_adjustment_decay_time        0:0:0
load_formula                      slots
schedd_job_info                   false
flush_submit_sec                  0
flush_finish_sec                  0
params                            none
reprioritize_interval             0:0:0
halftime                          168
usage_weight_list                 cpu=1.000000,mem=0.000000,io=0.000000
compensation_factor               5.000000
weight_user                       0.250000
weight_project                    0.250000
weight_department                 0.250000
weight_job                        0.250000
weight_tickets_functional         0
weight_tickets_share              0
share_override_tickets            TRUE
share_functional_shares           TRUE
max_functional_jobs_to_schedule   200
report_pjob_tickets               TRUE
max_pending_tasks_per_job         50
halflife_decay_list               none
policy_hierarchy                  OFS
weight_ticket                     0.010000
weight_waiting_time               0.000000
weight_deadline                   3600000.000000
weight_urgency                    0.100000
weight_priority                   1.000000
max_reservation                   0
default_duration                  INFINITY
SCHED
    our $pe_threaded = <<PE;
pe_name            threaded
slots              999
user_lists         NONE
xuser_lists        NONE
start_proc_args    NONE
stop_proc_args     NONE
allocation_rule    \$pe_slots
control_slaves     FALSE
job_is_first_task  TRUE
urgency_slots      min
accounting_summary TRUE
qsort_args         NONE
PE
    our $initialized = 0;
    
    sub periodic_method {
        return 'launch_extra_instances_and_terminate_old_instances';
    }
    
    sub on_exit_method {
        return 'terminate_all_instances';
    }
    
    around initialize_for_server {
        return if $initialized;
        $ec2_scheduler->initialize_for_server;
        
        mkdir($sge_confs_dir) unless -d $sge_confs_dir;
        
        # pe: give us a parallel environment for threaded jobs
        my $pe_file = file($sge_confs_dir, 'pe');
        unless (-s $pe_file) {
            open(my $fh, '>', $pe_file) || die "Could not write to $pe_file";
            print $fh $pe_threaded;
            close($fh);
        }
        system('qconf -Ap ' . $pe_file);
        
        # make sure we've configured the group and queue for all possible ec2
        # types
        foreach my $type (keys %VRPipe::Schedulers::ec2::instance_types) {
            # group
            my %host_groups;
            open(my $hgfh, 'qconf -shgrpl |') || die 'Could not open a pipe from qconf -shgrpl';
            while (<$hgfh>) {
                my ($group) = $_ =~ /^\@(\S+)/;
                $host_groups{$group} = 1 if $group;
            }
            if (exists $host_groups{allhosts}) {
                # get rid of the allhosts group
                system('qconf -dhgrp @allhosts');
            }
            unless (exists $host_groups{$type}) {
                # add a group for this type
                my $group_file = file($sge_confs_dir, 'hgrp.' . $type);
                unless (-s $group_file) {
                    open(my $fh, '>', $group_file) || die "Could not write to $group_file";
                    print $fh "group_name \@$type\nhostlist NONE\n";
                    close($fh);
                }
                system('qconf -Ahgrp ' . $group_file);
            }
            
            # queue
            my %queues;
            open(my $qfh, 'qconf -sql |') || die 'Could not open a pipe from qconf -sql';
            while (<$qfh>) {
                chomp;
                $queues{$_} = 1 if $_;
            }
            if (exists $queues{'all.q'}) {
                # get rid of the all.q queue
                system('qconf -dq all.q');
            }
            unless (exists $host_groups{$type}) {
                # add a queue for this type
                my $queue_file = file($sge_confs_dir, 'q.' . $type);
                unless (-s $queue_file) {
                    open(my $fh, '>', $queue_file) || die "Could not write to $queue_file";
                    my ($slots, $mb) = @{ $VRPipe::Schedulers::ec2::instance_types{$type} };
                    my $safe_mb = $mb - 100;
                    my $safe_gb = sprintf("%0.1f", $safe_mb / 1024);
                    my $safe_b  = $safe_mb * 1024 * 1024;
                    my $b       = $mb * 1024 * 1024;
                    print $fh "qname $type\nhostlist \@$type\n", $queue_defaults, "slots $slots\ncomplex_values slots=$slots,virtual_free=${safe_gb}G,mem_free=${safe_gb}G\ns_data $safe_b\nh_data $b\ns_rss $safe_b\nh_rss $b\ns_vmem $safe_b\nh_vmem $b\n";
                    close($fh);
                }
                system('qconf -Aq ' . $queue_file);
            }
        }
        
        # complex: enable memory and cpu reservation
        my $ca_file = file($sge_confs_dir, 'complex_attribs');
        unless (-s $ca_file) {
            open(my $fh, '>', $ca_file) || die "Could not write to $ca_file";
            print $fh $complex_attribs;
            close($fh);
        }
        system('qconf -Mc ' . $ca_file);
        
        # scheduler: make it "fill up hosts" instead of spreading jobs around
        # evenly, so that we're more likely to end up with empty nodes when
        # demand is low, so we can terminate them and save money
        my $s_file = file($sge_confs_dir, 'sconf');
        unless (-s $s_file) {
            open(my $fh, '>', $s_file) || die "Could not write to $s_file";
            print $fh $scheduling;
            close($fh);
        }
        system('qconf -Msconf ' . $s_file);
        
        # call sge's initialize_for_server method, if any
        $self->$orig();
        
        $initialized = 1;
    }
    
    method launch_extra_instances_and_terminate_old_instances (Str :$deployment!) {
        ## launch extra instances
        
        # look at what we have queuing to get a count of how many of each
        # instance type we need to launch
        open(my $qstatfh, 'qstat -g d -r -ne -s p |') || $self->throw('[sge_ec2scheduler] Unable to open a pipe from qstat -g d -r -ne -s p');
        my ($job_id, $slots, $array_index);
        my %types;
        my %job_id_to_type;
        my $minimum_queue_time = $deployment eq 'production' ? 300 : 30;
        
        # we're going to keep track of which job ids result in us launching new
        # instances, so that the next time this method gets called we won't
        # launch more instances for the same jobs (because SGE was slow at
        # starting them running on the new instances)
        
        while (<$qstatfh>) {
            #18 0.00000 vrpipe_174 ec2-user     qw    06/19/2013 14:23:01  2
            if (my ($this_job_id, $month, $day, $year, $hour, $min, $sec, $this_slots, $this_array_index) = $_ =~ /^\s+(\d+)\s+\S+\s+\S+\s+\S+\s+qw\s+(\d+)\/(\d+)\/(\d+)\s+(\d\d):(\d\d):(\d\d)\s+(\d+)(?:\s+(\d+))?/) {
                # we'll only consider jobs that have been pending for over 5mins
                my $dt = DateTime->new(
                    year      => $year,
                    month     => $month,
                    day       => $day,
                    hour      => $hour,
                    minute    => $min,
                    second    => $sec,
                    time_zone => 'local',
                );
                my $elapsed = time() - $dt->epoch;
                if ($elapsed >= $minimum_queue_time) {
                    $job_id      = $this_job_id;
                    $slots       = $this_slots;
                    $array_index = $this_array_index;
                }
                next;
            }
            
            if ($job_id) {
                $array_index ||= 0;
                my $lock_key = 'sge_ec2_job_id.' . $job_id . '.' . $array_index;
                next if $im->noted($lock_key);
                
                if (/^\s+Hard Resources:\s+\S+=(\d+)M/) {
                    my $mb = $1;
                    my $chosen_type;
                    my $typed = exists $job_id_to_type{$job_id};
                    if ($typed) {
                        $chosen_type = $job_id_to_type{$job_id};
                    }
                    else {
                        foreach my $type (@VRPipe::Schedulers::ec2::ordered_types) {
                            my ($available_slots, $available_mb) = @{ $VRPipe::Schedulers::ec2::instance_types{$type} };
                            $available_mb -= 500;
                            if ($slots <= $available_slots && $mb <= $available_mb) {
                                $chosen_type = $type;
                                last;
                            }
                        }
                    }
                    
                    if ($chosen_type) {
                        $types{$chosen_type} += $slots;
                        $im->note($lock_key, forget_after => 1800);
                        unless ($typed) {
                            $job_id_to_type{$job_id} = $chosen_type;
                        }
                    }
                    elsif (!$typed) {
                        $job_id_to_type{$job_id} = 0;
                    }
                }
            }
        }
        
        # launch the desired instance types
        my %launched_hosts;
        my $intitial_ec2_max_instances = $VRPipe::Schedulers::ec2::max_instances;
        while (my ($type, $count) = each %types) {
            # $count is the number of threads we want to run on this $type of
            # instance; adjust based on how many cores each instance has
            $count = ceil($count / $VRPipe::Schedulers::ec2::instance_types{$type}->[0]);
            $im->debug("[sge_ec2scheduler] Will launch $count $type instances");
            my @instances = $ec2_scheduler->launch_instances($type, $count);
            
            foreach my $instance (@instances) {
                my $pdn = $instance->privateDnsName;
                my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
                $launched_hosts{$host} = $type;
            }
        }
        
        if (keys %launched_hosts) {
            # start SGE on the new hosts
            my $tmp_config = $self->_altered_config('EXEC_HOST_LIST', join(' ', keys %launched_hosts));
            system('cd $SGE_ROOT; ./inst_sge -x -auto ' . $tmp_config) && $self->throw("[sge_ec2scheduler] Failed to run: inst_sge -x -auto $tmp_config");
            unlink($tmp_config);
            
            # add the new hosts to appropriate queues
            while (my ($host, $type) = each %launched_hosts) {
                system("qconf -aattr hostgroup hostlist $host \@$type") && $self->throw("[sge_ec2scheduler] Failed to run: qconf -aattr hostgroup hostlist $host \@$type");
            }
        }
        
        ## terminate old instances
        
        # start by seeing which hosts are not currently running any jobs
        open(my $qfh, 'qhost -j -ncb |') || $self->throw("[sge_ec2scheduler] Unable to open pipe from qhost -j -ncb");
        <$qfh>;
        <$qfh>;
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
        close($qfh) || $self->throw("[sge_ec2scheduler] Unable to close pipe from qhost -j -ncb");
        
        # consider empty hosts to be those that also haven't run a job in the
        # past 45mins
        my $min_uptime = 2700;
        my $max_do_nothing_time = $deployment eq 'production' ? $min_uptime : int($min_uptime / 10);
        my %to_terminate;
        my $own_pdn = $VRPipe::Schedulers::ec2::meta->privateDnsName;
        my ($own_host) = $own_pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
        while (my ($host, $details) = each %hosts) {
            # ones we just now launched probably won't be running anything yet
            next if exists $launched_hosts{$host};
            
            # obviously it's not empty if it is running a job right now
            next if $details->{jobs};
            
            # since we're charged by the hour, we won't terminate any nodes
            # unless they've been up for at least 45mins (in production)
            my ($ip) = $host =~ /(\d+-\d+-\d+-\d+)/;
            $ip =~ s/-/./g;
            next if $im->noted('starting_instance.' . $ip);
            my ($instance) = $ec2->describe_instances({ 'private-ip-address' => $ip });
            next unless $instance;
            if ($deployment eq 'production') {
                next unless $instance->up_time >= $min_uptime;
            }
            
            # don't terminate ourselves - the server that calls this method
            # won't have any handlers running on it
            my $pdn = $instance->privateDnsName;
            my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
            next if $host eq $own_host;
            
            # see if a job has run on this host in the past 45mins
            next if VRPipe::Job->search({ host => $host, end_time => { '>=' => DateTime->from_epoch(epoch => time() - $max_do_nothing_time) } });
            
            # don't terminate an instance that has a handler running on it right
            # now
            next if $ec2_scheduler->_handler_processes($host);
            
            # this host has been sitting there doing nothing for ages
            $to_terminate{$host} = $instance;
        }
        
        $self->_clean_up_hosts(\%to_terminate);
        
        # if the EC2 scheduler set its class variable $max_instances due to
        # a discovered quota, we'll return a method name and the current value
        # of that so that vrpipe-server can set it in globally, instead of the
        # change being lost because this method was called in a fork
        my $current_ec2_max_instances = $VRPipe::Schedulers::ec2::max_instances;
        if ($current_ec2_max_instances != $intitial_ec2_max_instances) {
            return ('_set_ec2_max_instances', $current_ec2_max_instances);
        }
    }
    
    method _set_ec2_max_instances (Int $mi) {
        $VRPipe::Schedulers::ec2::max_instances = $mi;
    }
    
    method _clean_up_hosts (HashRef $to_terminate) {
        while (my ($host, $instance) = each %$to_terminate) {
            my $type = $instance->instanceType;
            
            # try to remove the host from SGE a couple of times until it works
            my $worked = 0;
            for (1 .. 3) {
                $worked = 0;
                
                foreach my $cmd ("qconf -dattr hostgroup hostlist $host \@$type", "qconf -de $host", "qconf -dconf $host") {
                    $im->debug("[sge_ec2scheduler] About to try and do {$cmd}\n");
                    my $exit_code = system($cmd);
                    if ($exit_code) {
                        $im->log("[sge_ec2scheduler] Tried to run {$cmd} but it exited with code $exit_code");
                        last;
                    }
                    else {
                        $worked++;
                    }
                }
                
                last if $worked == 3;
                sleep(5);
            }
            
            unless ($worked == 3) {
                $im->log("[sge_ec2scheduler] Could not remove $host from SGE prior to instance termination");
            }
            
            $im->log("[sge_ec2scheduler] Will terminate instance $host");
            $instance->terminate;
        }
    }
    
    method terminate_all_instances {
        # (this is only run when the testing server goes down, so shouldn't
        #  affect the production server)
        my @all_instances = $ec2->describe_instances({
                'image-id'            => $VRPipe::Schedulers::ec2::ami,
                'availability-zone'   => $VRPipe::Schedulers::ec2::availability_zone,
                'instance-state-name' => 'running'
            }
        );
        
        # select which ones to terminate
        my $own_pdn = $VRPipe::Schedulers::ec2::meta->privateDnsName;
        my ($own_host) = $own_pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
        my %to_terminate;
        foreach my $instance (@all_instances) {
            # don't terminate ourselves - the server that calls this method
            # won't have any handlers running on it
            my $pdn = $instance->privateDnsName;
            my ($host) = $pdn =~ /(ip-\d+-\d+-\d+-\d+)/;
            next if $host eq $own_host;
            
            # don't terminate an instance that has a handler running on it right
            # now - possibly spawned by a production server
            next if $ec2_scheduler->_handler_processes($host);
            
            $to_terminate{$host} = $instance;
        }
        
        $self->_clean_up_hosts(\%to_terminate);
    }
    
    method _altered_config (Str $key, Str $val) {
        my $ofile = file($sge_confs_dir, 'config_tmp');
        open(my $fho, '>', $ofile)           || $self->throw("[sge_ec2scheduler] Could not write to $ofile");
        open(my $fhi, '<', $sge_config_file) || $self->throw("[sge_ec2scheduler] Could not read from $sge_config_file");
        my $found = 0;
        while (<$fhi>) {
            if (/^$key/) {
                print $fho qq[$key="$val"\n];
                $found = 1;
            }
            else {
                print $fho $_;
            }
        }
        close($fhi);
        
        unless ($found) {
            print $fho qq[$key="$val"\n];
        }
        close($fho);
        return $ofile;
    }
}

1;
