
=head1 NAME

VRPipe::FarmServer - track vrpipe-servers running against each farm

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::FarmServer extends VRPipe::Persistent {
    use DateTime;
    use POSIX qw(ceil);
    
    our $dead_time = 60;
    
    has 'farm' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'only_ours' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has 'last_claim_time' => (
        is          => 'rw',
        isa         => Datetime,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'heartbeat' => (
        is      => 'rw',
        isa     => Datetime,
        coerce  => 1,
        traits  => ['VRPipe::Persistent::Attributes'],
        default => sub { DateTime->now() }
    );
    
    __PACKAGE__->make_persistent();
    
    method time_since_heartbeat {
        my $heartbeat = $self->heartbeat;
        my $t         = time();
        return $t - $heartbeat->epoch;
    }
    
    method alive {
        my $elapsed = $self->time_since_heartbeat;
        my $alive   = $elapsed < $dead_time;
        unless ($alive) {
            $self->delete;     # (see note in DEMOLISH)
        }
        return $alive;
    }
    
    method still_alive {
        $self->heartbeat(DateTime->now());
    }
    
    method claim_setups {
        my $farm = $self->farm;
        
        # delete any farms no longer alive
        $self->search_rs({ heartbeat => { '<' => DateTime->from_epoch(epoch => time() - $dead_time) } })->delete;
        
        my @setups_to_trigger;
        my $transaction = sub {
            # unless we're only doing setups assigned directly to us, we want to
            # share setups evenly amongst all the farms
            unless ($self->only_ours) {
                my @farms = $self->search({ only_ours => 0 }, { for => 'update', order_by => { -asc => 'last_claim_time' } });
                
                my @search = ({ active => 1, desired_farm => undef, controlling_farm => undef }, { for => 'update' });
                my $num_setups = VRPipe::PipelineSetup->search(@search);
                
                if ($num_setups) {
                    my $pager = VRPipe::PipelineSetup->search_paged(@search, ceil($num_setups / scalar(@farms)));
                    
                    foreach my $fs (@farms) {
                        my $setups = $pager->next;
                        if ($fs->farm eq $farm) {
                            foreach my $setup (@$setups) {
                                $setup->controlling_farm($farm);
                                $setup->update;
                                push(@setups_to_trigger, $setup);
                            }
                            last;
                        }
                    }
                    
                    # the usual case is that $num_setups will be 1, and we don't
                    # want to always have it be claimed by the first farm in the
                    # table, which is why we sorted by last_claim_time, and now
                    # need to actually set that time for us so we're least likely
                    # to get the claim the next time a single new setup turns up
                    $self->last_claim_time(DateTime->now());
                    $self->update;
                }
            }
            
            # obviously we also claim setups assigned to us
            foreach my $setup (VRPipe::PipelineSetup->search({ active => 1, desired_farm => $farm, controlling_farm => undef }, { for => 'update' })) {
                $setup->controlling_farm($farm);
                $setup->update;
                push(@setups_to_trigger, $setup);
            }
        };
        $self->do_transaction($transaction, "Could not claim PipelineSetups for farm $farm");
        
        # since @setups_to_trigger were previously uncontrolled setups they are
        # most likely brand new and have never been run before, which means
        # there are no submissions for the farm server to manage; we'll start
        # off proceedings with a full trigger on all dataelements
        foreach my $setup (@setups_to_trigger) {
            $setup->trigger();
        }
        
        return VRPipe::PipelineSetup->search({ active => 1, controlling_farm => $farm });
    }
    
    sub DEMOLISH {
        my $self = shift;
        $self->delete;
        # we could do something like remove our farm name from the
        # controlling_farm of all PipelineSetups, but probably bad things would
        # happen if a setup that was still running was suddenly controlled by
        # a different farm. So once claimed, a setup will always run only on
        # that farm, even when its server goes down and doesn't come back
    }
}

1;
