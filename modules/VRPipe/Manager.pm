
=head1 NAME

VRPipe::Manager - methods for managing the execution of pipelines

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This is the main module used to discover what work needs to be done ('trigger'
each L<VRPipe::PipelineSetup>) and then to 'dispatch' that work out to the
system's job scheduler to actually run the command lines on the compute
cluster.

The B<vrpipe-trigger_pipelines> and B<vrpipe-dispatch_pipelines> scripts call
methods of this module.

Note that the current implementation is slow and inefficient, with lots of
serial looping over thousands of objects. A radical overhaul is planned.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Manager extends VRPipe::Persistent {
    use Parallel::ForkManager;
    use Sys::CPU;
    use POSIX qw(ceil);
    use Sys::Hostname::Long;
    
    our $DEFAULT_MAX_PROCESSES = Sys::CPU::cpu_count();
    our %step_limits;
    our %setups_with_step_limits;
    our %setups_checked_for_step_limits;
    our %step_limit_instigators;
    our $do_step_limits = 0;
    
    has 'global_limit' => (
        is  => 'rw',
        isa => PositiveInt
    );
    
    __PACKAGE__->make_persistent();
    
    around get (ClassName|Object $self: Persistent :$id?, PositiveInt :$global_limit = 500) {
        return $self->create(global_limit => $global_limit);
    }
    
    around create (ClassName|Object $self: Persistent :$id?, PositiveInt :$global_limit = 500) {
        # our db-storage needs are class-wide, so we only have 1 row in our
        # table
        return $self->$orig(id => 1, global_limit => $global_limit);
    }
    
    method register_farm_server (Str $farm, Bool :$only_ours = 0, Bool :$die_when_murdered = 0) {
        my $transaction = sub {
            my ($fs) = VRPipe::FarmServer->search({ farm => $farm }, { for => 'update' });
            return if ($fs && $fs->alive);
            $fs = VRPipe::FarmServer->create(farm => $farm, hostname => hostname_long, only_ours => $only_ours, die_when_murdered => $die_when_murdered);
            $fs->start_beating;
            return $fs;
        };
        
        return $self->do_transaction($transaction, "Could not register farm $farm");
    }
    
    method setups {
        return VRPipe::PipelineSetup->search({});
    }
}

1;
