=head1 NAME

VRPipe::StepAdaptor - connect the output of a step to the input of another

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A StepAdaptor defines which output of a step (or the DataElement) in a
L<VRPipe::Pipeline> should be passed to which other L<VRPipe::Step>(s) in the
same Pipeline.

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

class VRPipe::StepAdaptor extends VRPipe::Persistent {
    has 'pipeline' => (is => 'rw',
                       isa => Persistent,
                       coerce => 1,
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       belongs_to => 'VRPipe::Pipeline');
    
    has 'to_step' => (is => 'rw',
                      isa => IntSQL[4],
                      traits => ['VRPipe::Persistent::Attributes'],
                      is_key => 1);
    
    has 'adaptor_hash' => (is => 'rw',
                           isa => 'HashRef',
                           traits => ['VRPipe::Persistent::Attributes'],
                           default => sub { {} });
    
    __PACKAGE__->make_persistent();
    
    method adapt (Str :$input_key!, PreviousStepOutput :$pso?, VRPipe::DataElement :$data_element?) {
        if ((defined $pso ? 1 : 0) + (defined $data_element ? 1 : 0) != 1) {
            $self->throw("Exactly one of pso or data_element must be supplied");
        }
        
        my $hash = $self->adaptor_hash;
        
        if (defined $hash->{$input_key}) {
            my $ref = $hash->{$input_key};
            
            my @results;
            while (my ($key, $from_step) = each %$ref) {
                if ($data_element) {
                    if ($key eq 'data_element' && $from_step == 0) {
                        my $result = $data_element->result;
                        my $paths = $result->{paths} || $self->throw("data element ".$data_element->id." gave a result with no paths");
                        foreach my $path (@$paths) {
                            push(@results, VRPipe::File->get(path => file($path)->absolute));
                        }
                    }
                }
                else {
                    my $result = $pso->{$key}->{$from_step} || next;
                    push(@results, (ref($result) && ref($result) eq 'ARRAY') ? @$result : $result);
                }
            }
            
            return @results ? \@results : undef;
        }
        
        return;
    }
}

1;