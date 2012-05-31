=head1 NAME

VRPipe::StepOption - a user option for a Step

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

When defining a L<VRPipe::Step>, associating one or more of these will present
the user (when they are setting up their Pipeline that uses the Step) with the
option. In the C<body_sub()> of the Step the user's chosen options are available
with the C<options()> method from L<VRPipe::StepRole>.

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

class VRPipe::StepOption extends VRPipe::Persistent {
    has 'description' => (is => 'rw',
                          isa => Varchar[255],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    has 'optional' => (is => 'rw',
                       isa => 'Bool',
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       default => 0,
                       allow_key_to_default => 1);
    
    has 'default_value' => (is => 'rw',
                            isa => Varchar[255],
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            default => '',
                            allow_key_to_default => 1);
    
    has 'allowed_values' => (is => 'rw',
                            isa => 'ArrayRef',
                            traits => ['VRPipe::Persistent::Attributes'],
                            is_key => 1,
                            default => sub { [] },
                            allow_key_to_default => 1);
    
    __PACKAGE__->make_persistent();
}

1;