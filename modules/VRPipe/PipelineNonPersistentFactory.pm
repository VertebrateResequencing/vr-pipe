
=head1 NAME

VRPipe::PipelineNonPersistentFactory - a factory for pipelines defined in files

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

L<VRPipe::Pipeline>s are objects stored in the B<VRPipe> database. But
pipelines can also be defined in Perl module files. These non-persistent
versions are read from the .pm files and converted into Persistent objects as
necessary. In this way we can can create .pm files and also create pipelines on
the fly using the Persistent API, and both seem to behave the same way.

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

package VRPipe::PipelineNonPersistentFactory;
use MooseX::AbstractFactory;

implementation_does qw/VRPipe::PipelineRole/;
implementation_class_via sub { 'VRPipe::Pipelines::' . shift };

1;
