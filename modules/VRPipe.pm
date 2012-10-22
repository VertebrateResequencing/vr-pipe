
=head1 NAME

VRPipe - use VRPipe modules in your own Perl scripts

=head1 SYNOPSIS
<<<<<<< Updated upstream
    
    perl -MVRPipe -e '$vrfile = VRPipe::File->get(path => shift); ...' /tmp/foo
    perl -MVRPipe=testing -e '...'
    
=======
    
    perl -MVRPipe -e '$vrfile = VRPipe::File->get(path => shift); ...' /tmp/foo
    perl -MVRPipe=testing -e '...'

>>>>>>> Stashed changes     use VRPipe 'production'; # same as 'use VRPipe;'   
 use VRPipe 'testing';

=head1 DESCRIPTION

To use VRPipe in your scripts you have to choose the desired deployment and
load the schema. This module provides a convenient way of doing this. Simply
<<<<<<< Updated upstream
using it loads all VRPipe::Persistent modules (which account for the majority
of VRPipe modules), so you probably won't have to 'use' any others. Using the
"import" syntax you can also easily choose your deployment (production, which
is the default, or testing).
=======
using it loads all VRPipe::Persistent modules (which account for the majority of
VRPipe modules), so you probably won't have to 'use' any others. Using the
"import" syntax you can also easily choose your deployment (production, which is
the default, or testing).
>>>>>>> Stashed changes

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

class VRPipe extends VRPipe::Persistent::Schema {
    sub import {
        my ($self, $deployment) = @_;
        $deployment ||= 'production';
        unless ($deployment eq 'production' || $deployment eq 'testing') {
            $self->throw("'$deployment' is an invalid deployment: must be one of 'production' or 'testing'");
        }
        
        VRPipe::Persistent::SchemaBase->database_deployment($deployment);
    }
}

1;
