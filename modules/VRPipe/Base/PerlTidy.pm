
=head1 COPYRIGHT AND LICENSE

Author: Sendu Bala <sb10@sanger.ac.uk>.

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

package VRPipe::Base::PerlTidy;
use Perl::Tidy qw();
use strict;
use warnings;

sub perltidy {
    return Perl::Tidy::perltidy(prefilter => \&prefilter, postfilter => \&postfilter, @_);
}

sub prefilter {
    $_ = $_[0];
    
    # turn method into sub
    s/^(\s*)method (.*)/$1sub $2 \#__METHOD/gm;
    
    # turn class into simple braced block
    s/^(\s*)class (.+?)\{(.*?)\n/$1\{ \#__CLASS $2 \#__EXTRA$3\n/gm;
    
    # it messes up indentation and syntax checking for add_method(); turn it
    # into a plain sub
    s/\n(\s*)\$meta->add_method\('(\w+?)' => sub \{(.*?)\n(.+?)\n(\s*)\}\);\n/\n${1}sub $2 \{ \#__ADDMETHOD $3\n$4\n$5\}\n/gs;
    
    # it fails with syntax error on try/catch blocks, and can't get the same
    # solution as method blocks to work right; turn them into if/else
    s/\n(\s*)try \{(.+)?\}(\s+)catch ([^\{]*)\{/\n${1}if (1) { #__TRY $2}${3}else { #__CATCH $4 __ENDCATCH/gs;
    
    return $_;
}

sub postfilter {
    $_ = $_[0];
    
    # turn sub back into method
    s/^(\s*)sub (.*?)\s* \#__METHOD/${1}method $2/gm;
    
    # restore class
    s/^(\s*)\{\s+\#__CLASS (.+?) \#__EXTRA(.*?)\n/${1}class $2\{$3\n/gm;
    
    # restore add_method
    s/\n(\s*)sub (\w+?) \{ +\#__ADDMETHOD(.*?)\n(.+?)\n(\s*)\}\n/\n${1}\$meta->add_method\('$2' => sub \{ $3\n$4\n$5\}\);\n/gs;
    
    # restore try/catch
    s/\n(\s*)if \(1\) \{ +\#__TRY(.+)?\}(\s+)else \{ +\#__CATCH ([^\n]*) __ENDCATCH\n/\n${1}try \{$2\}${3}catch $4\{\n/gs;
    
    # the method->sub->method trick screws up comments that appear on the next
    # line; fix them now
    s/^(\s*method .+)\n\s+(\#.+?)\n(\s*)/$1\n$3$2\n$3/gm;
    
    # it also likes adding extra spaces between ; and # - remove these
    s/^(\s*\S[^#]+?)([;,]) +\#/$1$2 #/gm;
    
    # don't have completely empty lines after start of brace block
    s/\{\n\n/\{\n/g;
    
    # don't have completely empty lines between things; instead indent to the
    # correct level
    s/\n\n([\s]+)(\S.*)/\n$1\n$1$2/gm;
    
    return $_;
}

1;
