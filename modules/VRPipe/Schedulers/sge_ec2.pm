
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

Note that it takes alters configuration of your SGE installation, changing
things like your complex attributes, queues and groups.

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
    
    my $vrp_config = VRPipe::Config->new();
    our $sge_config_file = $vrp_config->sge_config_file;
    my $log_dir_method = VRPipe::Persistent::SchemaBase->database_deployment . '_logging_directory';
    my $sge_confs_dir = dir($vrp_config->$log_dir_method(), '.sge_confs');
    our $ec2_scheduler = VRPipe::SchedulerMethodsFactory->create('ec2', {});
    no warnings 'once';
    our $ec2 = $VRPipe::Schedulers::ec2::ec2;
    
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
pe_list               make
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
    
    # make sure we've configured the group and queue for all possible ec2
    # types
    mkdir($sge_confs_dir) unless -d $sge_confs_dir;
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
                my $safe_mb = $mb - 500;
                my $safe_gb = sprintf("%0.1f", $safe_mb / 1024);
                my $safe_b  = $safe_mb * 1024 * 1024;
                my $b       = $mb * 1024 * 1024;
                print $fh "qname $type\nhostlist \@$type\n", $queue_defaults, "slots $slots\ncomplex_values slots=$slots,virtual_free=${safe_gb}G,mem_free=${safe_gb}G\ns_data $safe_b\nh_data $b\ns_rss $safe_b\nh_rss $b\ns_vmem $safe_b\nh_vmem $b\n";
                close($fh);
            }
            system('qconf -Aq ' . $queue_file);
        }
    }
    
    # complex
    my $ca_file = file($sge_confs_dir, 'complex_attribs');
    unless (-s $ca_file) {
        open(my $fh, '>', $ca_file) || die "Could not write to $ca_file";
        print $fh $complex_attribs;
        close($fh);
    }
    system('qconf -Mc ' . $ca_file);
    
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
        my %launched_hosts;
        
        if (keys %launched_hosts) {
            # start SGE on the new hosts
            my $tmp_config = $self->_altered_config('EXEC_HOST_LIST', join(' ', keys %launched_hosts));
            system('cd $SGE_ROOT; ./inst_sge -x -auto ' . $tmp_config) && $self->throw("Failed to run: inst_sge -x -auto $tmp_config");
            unlink($tmp_config);
            
            # add the new hosts to appropriate queues
            while (my ($host, $type) = each %launched_hosts) {
                system('qconf -aattr hostgroup hostlist $host \@$type') && $self->throw("Failed to run: qconf -aattr hostgroup hostlist $host \@$type");
            }
        }
        
        ## terminate old instances
        
        # start by seeing which hosts are not currently running any jobs
        open(my $qfh, 'qhost -j -ncb |') || $self->throw("Unable to open pipe from qhost -j -ncb");
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
        
        if (@empty_hosts) {
            # stop SGE on the empty hosts
            my $tmp_config = $self->_altered_config('EXEC_HOST_LIST_RM', join(' ', @empty_hosts));
            #[edit /gluster/vrpipe/.sge_config to set host in EXEC_HOST_LIST_RM]
            system('cd $SGE_ROOT; ./inst_sge -ux -auto ' . $tmp_config) && $self->throw("Failed to run: inst_sge -ux -auto $tmp_config");
            unlink($tmp_config);
            
            # terminate the instances
            my $dt_parser = DateTime::Format::Natural->new;
            my $own_pdn   = $VRPipe::Schedulers::ec2::meta->privateDnsName;
            foreach my $host (@empty_hosts) {
                my ($ip) = $host =~ /(\d+)-(\d+)-(\d+)-(\d+)/;
                $ip =~ s/-/./;
                my ($instance) = $ec2->describe_instances({ 'private-ip-address' => $ip });
                
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
    }
    
    method terminate_all_instances {
        $ec2_scheduler->terminate_all_instances;
    }
    
    method _altered_config (Str $key, Str $val) {
        my $ofile = file($sge_confs_dir, 'config_tmp');
        open(my $fho, '>', $ofile)           || $self->throw("Could not write to $ofile");
        open(my $fhi, '<', $sge_config_file) || $self->throw("Could not read from $sge_config_file");
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
