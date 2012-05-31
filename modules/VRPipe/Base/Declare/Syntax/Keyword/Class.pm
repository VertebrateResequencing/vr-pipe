=head1 NAME

VRPipe::Base::Declare::Syntax::Keyword::Class - modify MooseX::Declare

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

use MooseX::Declare;

class VRPipe::Base::Declare::Syntax::Keyword::Class extends MooseX::Declare::Syntax::Keyword::Class {
    around import_symbols_from (Object $ctx) {
        'VRPipe::Base::UseMoose';
    }
    
    #before add_namespace_customizations (Object $ctx, Str $package) {
    #    $ctx->add_preamble_code_parts(
    #        'use metaclass (metaclass => "Moose::Meta::Class");' #, error_class => "MooseX::Error::Exception::Class"  *** MooseX::Error::Exception::Class does not actually compile
    #    );
    #}
    
    after add_namespace_customizations (Object $ctx, Str $package) {
        $ctx->add_preamble_code_parts(
            'use MooseX::StrictConstructor; use VRPipe::Base::Types ":all"; use TryCatch; use Path::Class;'
        );
    }
}

1;