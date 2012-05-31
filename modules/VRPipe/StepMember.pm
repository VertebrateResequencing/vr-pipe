=head1 NAME

VRPipe::StepMember - a Step of a Pipeline

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A L<VRPipe::Step> can be reused in multiple different L<VRPipe::Pipeline>s.
A Pipeline itself does not hold a list of Steps within its definition. Instead
it C<has_many> StepMembers, which tie a Step to a particular Pipeline and also
notes the C<step_number()>, so that the correct order of Steps is known.

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

class VRPipe::StepMember extends VRPipe::Persistent {
    has 'step' => (is => 'rw',
                   isa => Persistent,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1,
                   belongs_to => 'VRPipe::Step');
    
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'step_number' => (is => 'rw',
                          isa => IntSQL[4],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    __PACKAGE__->make_persistent();
    
    around step (PreviousStepOutput :$previous_step_outputs?, VRPipe::StepState :$step_state?) {
        my $step = $self->$orig();
        $step->step_state($step_state) if $step_state;
        $step->previous_step_outputs($previous_step_outputs) if $previous_step_outputs;
        return $step;
    }
}

1;