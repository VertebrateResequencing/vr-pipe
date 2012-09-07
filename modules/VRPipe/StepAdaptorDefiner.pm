
=head1 NAME

VRPipe::StepAdaptorDefiner - a non-persistent definer of StepAdaptors

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

L<VRPipe::StepAdaptor>s are Persistent objects stored in the database, and are
also a bit difficult to specify in code. This class makes it easy to define one
in a C<VRPipe::Pipeline::[pipeline_name]> module file, and it will
automatically be made into a real StepAdaptor.

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

class VRPipe::StepAdaptorDefiner {
    has 'from_step' => (
        is  => 'ro',
        isa => 'Int'
    );
    
    has 'to_step' => (
        is     => 'ro',
        isa    => PositiveInt,
        coerce => 1
    );
    
    has 'from_key' => (
        is      => 'ro',
        isa     => 'Str',
        builder => '_from_key_builder',
        lazy    => 1
    );
    
    has 'to_key' => (
        is  => 'ro',
        isa => 'Str'
    );
    
    method _from_key_builder {
        if ($self->from_step == 0) {
            return 'data_element';
        }
        else {
            $self->throw("from_key must be supplied if from_step is not 0");
        }
    }
    
    method define (Persistent|VRPipe::Pipeline $pipeline) {
        my $sa = VRPipe::StepAdaptor->create(pipeline => $pipeline, to_step => $self->to_step);
        my $adaptor_hash = $sa->adaptor_hash;
        $adaptor_hash->{ $self->to_key }->{ $self->from_key } = $self->from_step;
        $sa->adaptor_hash($adaptor_hash);
        $sa->update;
        return $sa;
    }
}

1;
