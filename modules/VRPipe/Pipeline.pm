
=head1 NAME

VRPipe::Pipeline - describes a pipeline (a series of steps)

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A Pipeline C<has_many> L<VRPipe::StepMember>s. It also has a name and
description. It also has associated objects that define how data flows in and
out of its component steps, and has behaviours for what should happen after
steps finish. Pipeline has a key place in the object model, being one of the
main objects users are likely to think about. Users want to "run a pipeline" on
their data. The pipelines that they can run are defined by these objects. A
Pipeline doesn't actually "do" anything, it just defines what is supposed to
happen and in what order.

It is possible to create a pipeline on the fly by C<get()>ing one of these
objects (suppling a new unique name), and then associating existing or new
StepMembers, StepAdaptors and StepBehvaiours. This is how the
B<vrpipe-pipeline_create> frontend lets users create their own pipelines.

Normally, however, you would define a Pipeline in a .pm file as
C<VRPipe::Pipelines::[your_pipeline_name]>. Special tricks automatically
convert these into VRPipe::Pipeline (and associated) objects in the database.

You must be wary, however, that the tricks cannot cope with the order of steps
changing, steps being removed or inserted. Only adding an extra step on the end
of a pipeline is safe. The advice is, test and stabalise a Pipeline, and once
it is in production never change it. If it turns out there are changes to steps
that have to be made, stop using that pipeline and create a new one with a
different name.

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

class VRPipe::Pipeline extends VRPipe::Persistent with VRPipe::PipelineRole {
    has 'name' => (is     => 'rw',
                   isa    => Varchar [64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has '_num_steps' => (is      => 'rw',
                         isa     => IntSQL [4],
                         traits  => ['VRPipe::Persistent::Attributes'],
                         default => 0);
    
    has 'description' => (is          => 'rw',
                          isa         => Text,
                          traits      => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    __PACKAGE__->make_persistent(has_many => [steps => 'VRPipe::StepMember']);
    
    # steps must be called to initially create stepmembers for a new pipeline,
    # but currently causes a memory leak, so we want to call it only once.
    # Everywhere else we call step_members to just get back the existing
    # step members
    sub step_members {
        my $self = shift;
        return VRPipe::StepMember->search({ pipeline => $self->id });
    }
}

1;
