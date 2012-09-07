
=head1 NAME

VRPipe::Base::Configuration::Env - allows configure to use environement vars

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

class VRPipe::Base::Configuration::Env {
    use overload(q[""] => 'stringify', fallback => 1);
    __PACKAGE__->meta->add_package_symbol('&()'  => sub { });
    __PACKAGE__->meta->add_package_symbol('&(""' => sub { shift->stringify });
    
    has variable => (
        is       => 'ro',
        isa      => 'Str',
        required => 1
    );
    
    # we don't have a value attribute so that the value is never stored in the
    # config file
    
    method value {
        return $ENV{ $self->variable };
    }
    
    method stringify {
        $self->value ? $self->value : '';
    }
}
