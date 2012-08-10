
=head1 NAME

VRPipe::LocalSchedulerJob - job tracking for the local scheduler

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

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

class VRPipe::LocalSchedulerJob extends VRPipe::Persistent {
    has 'cmd' => (is     => 'rw',
                  isa    => Text,
                  traits => ['VRPipe::Persistent::Attributes']);
    
    has 'cwd' => (is     => 'rw',
                  isa    => Dir,
                  coerce => 1,
                  traits => ['VRPipe::Persistent::Attributes']);
    
    has 'env' => (is     => 'rw',
                  isa    => 'HashRef',
                  traits => ['VRPipe::Persistent::Attributes']);
    
    has 'array_size' => (is      => 'rw',
                         isa     => IntSQL [8],
                         traits  => ['VRPipe::Persistent::Attributes'],
                         default => 1);
    
    has 'creation_time' => (is      => 'rw',
                            isa     => Datetime,
                            coerce  => 1,
                            traits  => ['VRPipe::Persistent::Attributes'],
                            default => sub { DateTime->now() });
    
    __PACKAGE__->make_persistent(has_many => [jobstates => 'VRPipe::LocalSchedulerJobState']);
}
