
=head1 NAME

VRPipe::Steps::r - a step

=head1 DESCRIPTION

Generic step for steps using R command, providing params for the R executable
and R libaries pathset, and generating a standard R command prefix

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::r extends VRPipe::Steps::r_script {
    has 'r_bin_path' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    method _build_standard_options {
        return ['r_bin_path', 'r_libs'];
    }
    
    method r_cmd_prefix {
        return $self->_cmd_prefix . $self->r_bin_path;
    }
    
    around options_definition {
        my %opts = %{ $self->$orig };
        delete $opts{rscript_cmd};
        return {
            %opts,
            r_bin_path => VRPipe::StepOption->create(description => 'path to your R executable and default arguments', optional => 1, default_value => 'R')
        };
    }
}

1;
