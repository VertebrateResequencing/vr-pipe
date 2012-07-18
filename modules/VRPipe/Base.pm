
=head1 NAME

VRPipe::Base - a set up module for all VRPipe classes

=head1 SYNOPSIS
    
    use VRPipe::Base;
    
    class VRPipe::MyClass {
        has 'my_attribute' => (is => 'rw', isa => VarChar[16]);
        
        method my_method (Str $input) {
            # throw an exception
            $self->throw("This is an exception") if $bad_thing_happened;
            
            # catch and handle a possible exception
            eval {
                my $other = VRPipe::othermodule->new()
                $other->dangerous_method();
            };
            if ($@) {
                # got an exception, do something about it...
            }
            else {
                # didn't get an exception, do something else...
            }
            
            # warn about something
            $self->warn("This is a warning") if $something_not_quite_right;
            
            # print a debugging message for those who want to read it
            $self->debug("This is a debugging message");
        }
    }

=head1 DESCRIPTION

Used in place of 'use MooseX::Declare;', providing the same functionality. In
addition it sets the base class of your class to VRPipe::Base:UseMoose, which
in turn gives you the VRPipe::Base::Debuggable and Cleanable roles.

Furthermore, all VRPipe::Base::Types are loaded for you, and
MooseX::StrictConstructor is turned on.

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

use MooseX::Declare;

class VRPipe::Base extends MooseX::Declare is dirty {
    use aliased 'VRPipe::Base::Declare::Syntax::Keyword::Class', 'ClassKeyword';
    use aliased 'VRPipe::Base::Declare::Syntax::Keyword::Role',  'RoleKeyword';
    use aliased 'MooseX::Declare::Syntax::Keyword::Namespace',   'NamespaceKeyword';
    
    clean;
    
    around keywords (ClassName $self:) {
        ClassKeyword->new(identifier => 'class'), RoleKeyword->new(identifier => 'role'), NamespaceKeyword->new(identifier => 'namespace'),;
    }
}

1;
