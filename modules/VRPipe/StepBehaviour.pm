
=head1 NAME

VRPipe::StepBehaviour - specify actions to take after a Pipeline Step completes

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

It is frequently useful to trigger some behaviour after a L<VRPipe::Step> in a
L<VRPipe::Pipeline> completes, but it wouldn't be appropriate to do this thing
as part of the Step itself, and you also want to make it optional with the user
deciding if it should happen. For example, Steps 2 and 3 of your Pipeline might
create files, both of which are needed by Step 4, but once Step 4 has completed
those files are usually no longer needed. You would want a behaviour that
triggered after Step 4 completes that deletes the output of Steps 2 and 3.

A StepBehaviour allows you to specify these kind of Pipeline-specific optional
behaviours. Behaviours are defined as methods in this class; more behaviours
should be added here.

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

class VRPipe::StepBehaviour extends VRPipe::Persistent {
    has 'pipeline' => (is         => 'rw',
                       isa        => Persistent,
                       coerce     => 1,
                       traits     => ['VRPipe::Persistent::Attributes'],
                       is_key     => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'after_step' => (is     => 'rw',
                         isa    => IntSQL [4],
                         traits => ['VRPipe::Persistent::Attributes'],
                         is_key => 1);
    
    has 'behaviour' => (is      => 'ro',
                        isa     => Varchar [64],
                        traits  => ['VRPipe::Persistent::Attributes'],
                        is_key  => 1,
                        default => 'delete_outputs');
    
    has 'behaviour_array' => (is      => 'rw',
                              isa     => 'ArrayRef',
                              traits  => ['VRPipe::Persistent::Attributes'],
                              default => sub { [] });
    
    has 'regulated_by' => (is      => 'rw',
                           isa     => Varchar [64],
                           traits  => ['VRPipe::Persistent::Attributes'],
                           default => '');
    
    has 'default_regulation' => (is      => 'rw',
                                 isa     => 'Bool',
                                 traits  => ['VRPipe::Persistent::Attributes'],
                                 default => 0);
    
    has 'dataelement' => (is  => 'rw',
                          isa => 'VRPipe::DataElement');
    has 'pipelinesetup' => (is  => 'rw',
                            isa => 'VRPipe::PipelineSetup');
    
    __PACKAGE__->make_persistent();
    
    method behave (VRPipe::DataElement :$data_element!, VRPipe::PipelineSetup :$pipeline_setup!) {
        my $go_ahead     = 0;
        my $regulated_by = $self->regulated_by;
        if ($regulated_by) {
            my $options = $pipeline_setup->options;
            if (defined $options->{$regulated_by}) {
                $go_ahead = $options->{$regulated_by};
            }
            else {
                $go_ahead = $self->default_regulation;
            }
        }
        else {
            $go_ahead = $self->default_regulation;
        }
        
        return unless $go_ahead;
        
        $self->dataelement($data_element);
        $self->pipelinesetup($pipeline_setup);
        
        my $array = $self->behaviour_array;
        
        foreach my $behaviour (@$array) {
            my ($method, @steps) = @$behaviour;
            $self->throw("'$method' is not a valid behaviour") unless $self->can($method);
            
            if ($method eq 'delete_inputs') {
                $self->throw('delete_inputs is a special behaviour that only operate on step 0, i.e. files from the datasource') unless (@steps == 1 && $steps[0] == 0);
                $self->$method;
            }
            else {
                return $self->$method($self->step_numbers_to_states(\@steps));
            }
        }
    }
    
    method step_numbers_to_states (ArrayRef[PositiveInt] $steps) {
        my $pipeline = $self->pipeline;
        my $de       = $self->dataelement;
        my $ps       = $self->pipelinesetup;
        
        my %step_nums = map { $_ => 1 } @$steps;
        
        my @states;
        foreach my $stepm ($pipeline->step_members) {
            if (exists $step_nums{ $stepm->step_number }) {
                push(@states, VRPipe::StepState->get(stepmember => $stepm, pipelinesetup => $ps, dataelement => $de));
            }
        }
        
        return \@states;
    }
    
    method delete_outputs (ArrayRef[VRPipe::StepState] $states) {
        foreach my $state (@$states) {
            $state->unlink_output_files;
        }
    }
    
    method start_over (ArrayRef[VRPipe::StepState] $states) {
        foreach my $state (@$states) {
            $state->start_over;
        }
    }
    
    method delete_inputs {
        my $data_element = $self->dataelement;
        my $result       = $data_element->result;
        my $paths        = $result->{paths} || $self->throw("data element " . $data_element->id . " gave a result with no paths");
        foreach my $path (@$paths) {
            VRPipe::File->get(path => file($path)->absolute)->unlink;
        }
    }
}

1;
