
=head1 NAME

VRPipe::StepBehaviourDefiner - non-Persistent definition of StepBehaviours

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A L<VRPipe::StepBehaviour> is a Persistent object stored in the database. This
class makes it possible to define one in a
C<VRPipe::Pipelines::[pipeline_name]> module file. It will automatically become
a real StepBehaviour.

It is necessary because a real StepBehaviour can't be created without the real
Pipeline already existing in the database, which obviously isn't the case for a
new Pipeline you're defining in .pm file.

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

class VRPipe::StepBehaviourDefiner {
    use Data::Compare;
    
    has 'after_step' => (is  => 'ro',
                         isa => PositiveInt);
    
    has 'behaviour' => (is  => 'ro',
                        isa => 'Str');
    
    has 'act_on_steps' => (is  => 'ro',
                           isa => 'ArrayRef');
    
    has 'regulated_by' => (is  => 'ro',
                           isa => 'Str');
    
    has 'default_regulation' => (is  => 'ro',
                                 isa => 'Bool');
    
    method define (Persistent|VRPipe::Pipeline $pipeline) {
        my $sb = VRPipe::StepBehaviour->create(pipeline => $pipeline, after_step => $self->after_step, behaviour => $self->behaviour);
        my $array = $sb->behaviour_array;
        
        my $behaviour    = $self->behaviour;
        my @steps        = sort { $a <=> $b } @{ $self->act_on_steps };
        my $already_have = 0;
        foreach my $existing (@$array) {
            my ($this_b, @these_steps) = @$existing;
            next unless $this_b eq $behaviour;
            
            if (Compare(\@steps, \@these_steps)) {
                $already_have = 1;
                last;
            }
        }
        
        unless ($already_have) {
            push(@$array, [$behaviour, @steps]);
        }
        
        $sb->behaviour_array($array);
        $sb->regulated_by($self->regulated_by);
        $sb->default_regulation($self->default_regulation);
        $sb->update;
        return $sb;
    }
}

1;
