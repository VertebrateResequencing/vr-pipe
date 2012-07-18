
=head1 NAME

VRPipe::PipelineSetup - set up a pipeline

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The main thing that users want to do is "run a pipeline" on a given set of
data. Pipelines are modelled by L<VRPipe::Pipline> objects, and the set of data
is modelled by a L<VRPipe::DataSource>. A PipelineSetup relates the two and
also stores the user configuration of both, resulting in a named object
defining what the user wants to happen.

So users "set up" pipelines and get back a PipelineSetup. They, and the
B<VRPipe> system itself, look to these PipelineSetups to see what is supposed
to be run and how, on what data.

Multiple PipelineSetups can run at once, even working on the same Pipelines and
DataSources (presumably with different configuration options), and the system
ensures that there are no problems with similar work being done by different
PipelineSetups overwriting each other's files.

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

class VRPipe::PipelineSetup extends VRPipe::Persistent {
    has 'name' => (is     => 'rw',
                   isa    => Varchar [64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'datasource' => (is         => 'rw',
                         isa        => Persistent,
                         coerce     => 1,
                         traits     => ['VRPipe::Persistent::Attributes'],
                         is_key     => 1,
                         belongs_to => 'VRPipe::DataSource');
    
    has 'pipeline' => (is         => 'rw',
                       isa        => Persistent,
                       coerce     => 1,
                       traits     => ['VRPipe::Persistent::Attributes'],
                       is_key     => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'output_root' => (is     => 'rw',
                          isa    => Dir,
                          coerce => 1,
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    has 'options' => (is                   => 'rw',
                      isa                  => 'HashRef',
                      traits               => ['VRPipe::Persistent::Attributes'],
                      default              => sub { {} },
                      allow_key_to_default => 1,
                      is_key               => 1);
    
    has 'description' => (is          => 'rw',
                          isa         => Text,
                          traits      => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    has 'active' => (is      => 'rw',
                     isa     => 'Bool',
                     traits  => ['VRPipe::Persistent::Attributes'],
                     default => 1);
    
    has 'user' => (is      => 'rw',
                   isa     => Varchar [64],
                   traits  => ['VRPipe::Persistent::Attributes'],
                   default => 'vrpipe');
    
    __PACKAGE__->make_persistent(has_many => [states => 'VRPipe::StepState']);
    
    # because lots of frontends get the dataelementstates for a particular
    # pipelinesetup
    sub _des_search_args {
        my $self = shift;
        return ({ pipelinesetup => $self->id, 'dataelement.withdrawn' => 0 }, { join => 'dataelement', prefetch => 'dataelement' });
    }
    
    method dataelementstates_pager {
        return VRPipe::DataElementState->search_paged($self->_des_search_args);
    }
    
    method dataelementstates {
        return VRPipe::DataElementState->search($self->_des_search_args);
    }
}

1;
