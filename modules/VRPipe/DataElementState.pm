
=head1 NAME

VRPipe::DataElementState - tracks the state of a DataElement for a given Setup

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A L<VRPipe::DataElement> can be the input to multiple different
L<VRPipe::PipelineSetup>s; this object keeps track of how many steps in the
pipeline a DataElement has gotten through in a particular PipelineSetup.

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

class VRPipe::DataElementState extends VRPipe::Persistent {
    has 'pipelinesetup' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::PipelineSetup'
    );
    
    has 'dataelement' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::DataElement'
    );
    
    has 'completed_steps' => (
        is      => 'rw',
        isa     => IntSQL [4],
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    __PACKAGE__->make_persistent();
    
    method start_from_scratch (ArrayRef[PositiveInt] $step_numbers?) {
        my $do_our_steps = 0;
        my %step_numbers;
        if ($step_numbers && @$step_numbers > 0) {
            %step_numbers = map { $_ => 1 } @$step_numbers;
        }
        else {
            # by default, we'll start_over all steps that produced output files
            # unique to our own stepstates
            %step_numbers = map { $_ => 1 } $self->our_step_numbers;
        }
        
        $self->pipelinesetup->log_event("DataElementState->start_from_scratch called for steps " . join(", ", sort { $a <=> $b } keys %step_numbers), dataelement => $self->dataelement->id, record_stack => 1);
        
        # get all the stepstates made for our dataelement and pipeline and
        # start_over() them
        foreach my $ss (VRPipe::StepState->search({ dataelement => $self->dataelement->id, pipelinesetup => $self->pipelinesetup->id }, { prefetch => 'stepmember' })) {
            next unless exists $step_numbers{ $ss->stepmember->step_number };
            $ss->pipelinesetup->log_event("Calling StepState->start_over as part of a DataElementState->start_from_scratch", stepstate => $ss->id, dataelement => $ss->dataelement->id);
            $ss->start_over();
            
            # ss->start_over deletes submissions for this stepstate - or for the
            # stepstate that this ss had same_submissions_as. We want to avoid
            # the situation where this ss is left with the same_submissions_as
            # another ss that has no submissions, and may be for a setup that
            # depends on us and so will never get any submissions
            if ($ss->same_submissions_as) {
                # delete the ss, so that when PipelineSetup->trigger() creates
                # it again, it will either get its own submissions, or have the
                # same_submissions_as something that actually has submissions
                $ss->pipelinesetup->log_event("Deleting the StepState because it has same_submissions_as, and we just deleted the Submissions of the master StepState", stepstate => $ss->id, dataelement => $ss->dataelement->id);
                $ss->delete;
            }
            else {
                # if other stepstates have the same_submissions_as this ss,
                # since we just deleted all our submissions they are now
                # invalid; delete them so they'll be recreated in a valid way by
                # PipelineSetup->trigger() at some point in the future
                #*** except that, what if the trigger never happens, because the
                #    other_ss setup is deactivated or complete? We won't see
                #    all its output files with vrpipe-output, and other strange
                #    things might happen... don't know how to solve this...
                #foreach my $other_ss (VRPipe::StepState->search({ same_submissions_as => $ss->id })) {
                #    eval { $other_ss->delete; };
                #}
            }
        }
        
        # each start_over() call will have set completed_steps(0) on us, but
        # we'll call it explicitly incase that ever changes
        $self->completed_steps(0);
        $self->update;
        
        # If this data element was used as the source of another dataelement, we want to also restart those dataelements
        my @children;
        foreach my $link (VRPipe::DataElementLink->search({ pipelinesetup => $self->pipelinesetup->id, parent => $self->dataelement->id }, { prefetch => 'child' })) {
            push @children, $link->child->element_states;
        }
        foreach my $child (@children) {
            $child->pipelinesetup->log_event("DataElementState->start_from_scratch for des " . $self->id . " called start_from_scratch on us because we are its child", dataelement => $child->dataelement->id);
            $child->start_from_scratch();
        }
        
        $self->pipelinesetup->log_event("DataElementState->start_from_scratch set completed_steps to 0 and will now return", dataelement => $self->dataelement->id);
    }
    
    method our_step_numbers {
        # get all our stepstates
        my @step_states;
        my %ss_ids;
        foreach my $ss (VRPipe::StepState->search({ dataelement => $self->dataelement->id, pipelinesetup => $self->pipelinesetup->id })) {
            push(@step_states, $ss);
            $ss_ids{ $ss->id } = 1;
        }
        
        my %step_nums;
        foreach my $ss (@step_states) {
            my @ofiles = $ss->_output_files;
            my $ours   = 1;
            # go through all the output files of our step state
            foreach my $sof (@ofiles) {
                next if $sof->output_key eq 'temp';
                
                # check that all other step states that output this same file
                # are also our own step states
                my $pager = VRPipe::StepOutputFile->search_paged({ file => $sof->file->id, 'me.id' => { '!=' => $sof->id }, output_key => { '!=' => 'temp' } }, { prefetch => 'stepstate' });
                PLOOP: while (my $other_sofs = $pager->next) {
                    foreach my $other_sof (@$other_sofs) {
                        my $other_ss = $other_sof->stepstate;
                        unless (exists $ss_ids{ $other_ss->id }) {
                            $ours = 0;
                            last PLOOP;
                        }
                    }
                }
                
                last unless $ours;
            }
            
            $ours || next;
            $step_nums{ $ss->stepmember->step_number } = 1;
        }
        
        my @step_nums = sort { $a <=> $b } keys %step_nums;
        return @step_nums;
    }
}

1;
