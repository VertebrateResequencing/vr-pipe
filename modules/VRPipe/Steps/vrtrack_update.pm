
=head1 NAME

VRPipe::Steps::vrtrack_update - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::vrtrack_update with VRPipe::StepRole {
    # eval these so that test suite can pass syntax check on this module when
    # VRTrack is not installed
    eval "use VRTrack::Factory;";
    
    method options_definition {
        return { vrtrack_db => VRPipe::StepOption->create(description => 'the name of your VRTrack database (other connection settings are taken from the standard VRTrack environment variables)') };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub { };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Empty shell step for basing other vrtrack-related steps on";
    }
    
    method max_simultaneous {
        return 75;
    }
    
    method get_vrtrack (ClassName|Object $self: Str :$db!, Str :$mode = 'rw') {
        return VRTrack::Factory->instantiate(database => $db, mode => $mode) || $self->throw("Could not connect to the database '$db'");
    }
}

1;
