
=head1 NAME

VRPipe::Requirements - describes what system resources are needed to run a Job

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

In order to execute a command line on a compute cluster, some job scheduling
system (like LSF) will be in place to select the node on which to run the job,
and to choose which job gets run first. For optimal efficiency the scheduler
should be told the likely amount of memory, time and other system resources the
job will need.

This Requirements class holds these scheduling hints. A L<VRPipe::Submission>
then ties a Requirements to a L<VRPipe::Job>.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Requirements extends VRPipe::Persistent {
    has 'memory' => (
        is     => 'rw',
        isa    => IntSQL [5],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'time' => (
        is     => 'rw',
        isa    => IntSQL [9],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'cpus' => (
        is                   => 'rw',
        isa                  => IntSQL [5],
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => 1,
        is_key               => 1,
        allow_key_to_default => 1
    );
    
    has 'tmp_space' => (
        is                   => 'rw',
        isa                  => IntSQL [5],
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => 0,
        is_key               => 1,
        allow_key_to_default => 1
    );
    
    has 'local_space' => (
        is                   => 'rw',
        isa                  => IntSQL [5],
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => 0,
        is_key               => 1,
        allow_key_to_default => 1
    );
    
    has 'custom' => (
        is                   => 'rw',
        isa                  => 'HashRef',
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => sub { {} },
        is_key               => 1,
        allow_key_to_default => 1
    );
    
    __PACKAGE__->make_persistent();
    
    # we used to store time as hours, but now as seconds
    around time (Maybe[PositiveInt] $time?) {
        if ($time) {
            if ($time < 60) {
                $time *= 60 * 60;
            }
            return $self->$orig($time);
        }
        
        my $time = $self->$orig;
        if ($time < 60) {
            $time *= 60 * 60;
        }
        return $time;
    }
}

1;
