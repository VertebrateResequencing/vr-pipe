
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
    has 'pipelinesetup' => (is         => 'rw',
                            isa        => Persistent,
                            coerce     => 1,
                            traits     => ['VRPipe::Persistent::Attributes'],
                            is_key     => 1,
                            belongs_to => 'VRPipe::PipelineSetup');
    
    has 'dataelement' => (is         => 'rw',
                          isa        => Persistent,
                          coerce     => 1,
                          traits     => ['VRPipe::Persistent::Attributes'],
                          is_key     => 1,
                          belongs_to => 'VRPipe::DataElement');
    
    has 'completed_steps' => (is      => 'rw',
                              isa     => IntSQL [4],
                              traits  => ['VRPipe::Persistent::Attributes'],
                              default => 0);
    
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
        
        # get all the stepstates made for our dataelement and pipeline and
        # start_over() them
        foreach my $ss (VRPipe::StepState->search({ dataelement => $self->dataelement->id, pipelinesetup => $self->pipelinesetup->id }, { prefetch => 'stepmember' })) {
            next unless exists $step_numbers{ $ss->stepmember->step_number };
            $ss->start_over();
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
            $child->start_from_scratch();
        }
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
                my $pager = VRPipe::StepOutputFile->search_paged({ file => $sof->file->id }, { prefetch => 'stepstate' });
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
