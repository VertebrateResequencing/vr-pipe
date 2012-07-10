=head1 NAME

VRPipe::DataElement - an element of a DataSource

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A DataElement respresents an element of a DataSource, and is the input to one
or more pipeline steps.

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

class VRPipe::DataElement extends VRPipe::Persistent {
    has 'datasource' => (is => 'rw',
                         isa => Persistent,
                         coerce => 1,
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1,
                         belongs_to => 'VRPipe::DataSource');
    
    has 'result' => (is => 'rw',
                     isa => 'HashRef',
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'withdrawn' => (is => 'rw',
                        isa => 'Bool',
                        traits => ['VRPipe::Persistent::Attributes'],
                        default => 0);
    
    __PACKAGE__->make_persistent(); # has_many => [element_states => 'VRPipe::DataElementState'] doesn't work because of ordering issues?
    
    method element_states {
        return VRPipe::DataElementState->search({ dataelement => $self->id }, { prefetch => [qw(dataelement pipelinesetup)] });
    }
    
    method start_from_scratch (VRPipe::PipelineSetup $setup, ArrayRef[PositiveInt] $step_numbers?) {
        VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $self)->start_from_scratch($step_numbers ? $step_numbers : ());
    }
}

1;