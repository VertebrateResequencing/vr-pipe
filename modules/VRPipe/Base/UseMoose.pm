
=head1 NAME

VRPipe::Base::UseMoose - getting VRPipe to use Moose

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

package VRPipe::Base::UseMoose;

use Moose ();
use Moose::Exporter;
use Moose::Util::MetaRole;

Moose::Exporter->setup_import_methods(also => 'Moose');

sub init_meta {
    shift;
    my %args = @_;
    
    Moose->init_meta(%args, base_class => 'VRPipe::Base::Moose');
    
    # for some reason this is required for something like 'after BUILDALL' to
    # work in VRPipe::Base::Debuggable
    #Moose::Util::MetaRole::apply_metaroles(
    #    for             => $args{for_class},
    #    class_metaroles => {
    #        class => => ['MooseX::NonMoose::Meta::Role::Class'],
    #        constructor => ['MooseX::NonMoose::Meta::Role::Constructor'],
    #    },
    #);
    
    # roles can be set like this, or by using 'with' in VRPipe::Base::Moose;
    # not sure which is best?
    #Moose::Util::MetaRole::apply_base_class_roles(
    #    for   => $args{for_class},
    #    roles => ['VRPipe::Base::Debuggable'],
    #);
    
    return $args{for_class}->meta();
}

no Moose;

1;
