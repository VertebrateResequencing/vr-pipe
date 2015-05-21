
=head1 NAME

VRPipe::DataSourceVRPipeRole - role for the vrpipe datasource

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This role just allows the vrpipe datasource to override the _create_elements
method.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

role VRPipe::DataSourceVRPipeRole with VRPipe::DataSourceFilterRole {
    around _create_elements (ArrayRef $e_args!) {
        # the all() method would be screwed up if anything overrides _create
        # elements() to change the element_args, so we check if it changed
        # and add new connections between element result and parent
        if (exists $self->{result_to_linkargs}) {
            my %current_res_to_str;
            foreach my $e_arg (@$e_args) {
                my $lookup = delete $e_arg->{lookup};
                $self->throw("the lookup key was missing from element args") unless $lookup;
                
                my $res_to_str = $self->_res_to_str($e_arg->{result});
                $current_res_to_str{$res_to_str} = 1;
                
                if ($res_to_str ne $lookup) {
                    $self->{result_to_linkargs}->{$res_to_str} = $self->{result_to_linkargs}->{$lookup};
                }
            }
            
            foreach my $res_to_str (keys %{ $self->{result_to_linkargs} }) {
                delete $self->{result_to_linkargs}->{$res_to_str} unless exists $current_res_to_str{$res_to_str};
            }
        }
        
        $self->$orig($e_args);
    }
    
    sub _res_to_str {
        my ($self, $res) = @_;
        my $args = { result => $res };
        VRPipe::DataElement->_convert_result($args);
        return $args->{filelist} . '|' . $args->{keyvallist};
    }
}

1;
