
=head1 NAME

VRPipe::Persistent::ConverterFactory - a factory for SQL converteres

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

package VRPipe::Persistent::ConverterFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::Persistent::ConverterRole/;
implementation_class_via sub { 'VRPipe::Persistent::Converter::' . shift };

1;
