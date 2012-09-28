
=head1 NAME

VRPipe::PipelineRole - a role required by all piplines

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

role VRPipe::PipelineRole {
    use VRPipe::StepAdaptorDefiner;
    use VRPipe::StepBehaviourDefiner;
    use Module::Find;
    
    requires 'name';
    requires 'description';
    
    has 'num_steps' => (
        is      => 'rw',
        isa     => 'Int',
        lazy    => 1,
        builder => '_build_num_steps'
    );
    
    our %pipeline_modules; # we need to delay finding our pipeline modules until the last moment, in case of lib changes
    our $found_modules = 0;
    
    method _build_num_steps {
        if ($self->can('step_members')) {
            return scalar $self->step_members;
        }
        else {
            return scalar @{ $self->_step_list->[0] };
        }
    }
    
    method steps {
        if ($self->can('step_members')) {
            return [map { $_->step } $self->step_members];
        }
        elsif ($self->can('step_names')) {
            return [map { VRPipe::Step->get(name => $_) } $self->step_names];
        }
        else {
            $self->throw("Pipeline " . $self->name . " has no steps!");
        }
    }
    
    method adaptors (VRPipe::Pipeline $persistent?) {
        if ($self->can('step_members')) {
            return [VRPipe::StepAdaptor->search({ pipeline => $self->id })];
        }
        elsif ($self->can('adaptor_definitions') && $persistent) {
            return [map { VRPipe::StepAdaptorDefiner->new(%{$_})->define($persistent) } $self->adaptor_definitions];
        }
        else {
            return [];
        }
    }
    
    method behaviours (VRPipe::Pipeline $persistent?) {
        if ($self->can('step_members')) {
            return [VRPipe::StepBehaviour->search({ pipeline => $self->id })];
        }
        elsif ($self->can('behaviour_definitions') && $persistent) {
            return [map { VRPipe::StepBehaviourDefiner->new(%{$_})->define($persistent) } $self->behaviour_definitions];
        }
        else {
            return [];
        }
    }
    
    method add_step (VRPipe::Step $step) {
        my $step_num = $self->num_steps() + 1;
        my $sm = VRPipe::StepMember->create(step => $step, pipeline => $self, step_number => $step_num);
        $self->num_steps($step_num);
        return $sm;
    }
    
    method _construct_pipeline (ArrayRef[VRPipe::Step] $steps, ArrayRef[VRPipe::StepAdaptor] $adaptors, ArrayRef[VRPipe::StepBehaviour] $behaviours) {
        # first check the pipeline hasn't already been constructed correctly
        my $all_ok   = 1;
        my $step_num = 0;
        my @sms;
        foreach my $step (@$steps) {
            my $step_id = $step->id;
            my @results = VRPipe::StepMember->search({
                    'step'        => $step_id,
                    'pipeline'    => $self->id,
                    'step_number' => ++$step_num
                }
            );
            
            if (@results == 1) {
                push(@sms, @results);
            }
            else {
                $all_ok = 0;
            }
        }
        
        unless ($all_ok) {
            # delete all stepmembers currently associated with the pipeline, as
            # long as all stepstates associated with those stepmembers are
            # complete (else throw)
            # *** if you alter a pipeline mid-run, are you forced to restart from scratch?!
            #*** not yet implemented
            
            # construct pipeline from scratch
            $self->num_steps(0);
            @sms = ();
            foreach my $step (@$steps) {
                push(@sms, $self->add_step($step));
            }
        }
        
        # create adaptors, deleting ones we no longer want
        my %wanted_sas;
        foreach my $sa (@$adaptors) {
            $wanted_sas{ $sa->id } = 1;
        }
        foreach my $sa (VRPipe::StepAdaptor->search({ pipeline => $self->id })) {
            unless (exists $wanted_sas{ $sa->id }) {
                $sa->delete;
            }
        }
        
        # create behaviours, deleting ones we no longer want
        my %wanted_sbs;
        foreach my $sb (@$behaviours) {
            $wanted_sbs{ $sb->id } = 1;
        }
        foreach my $sb (VRPipe::StepBehaviour->search({ pipeline => $self->id })) {
            unless (exists $wanted_sbs{ $sb->id }) {
                $sb->delete;
            }
        }
        
        return @sms;
    }
}

1;
