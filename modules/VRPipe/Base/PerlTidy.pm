package VRPipe::Base::PerlTidy;
use Perl::Tidy qw();
use strict;
use warnings;

sub perltidy {
    return Perl::Tidy::perltidy( prefilter => \&prefilter, postfilter => \&postfilter, @_ );
}

sub prefilter {
    $_ = $_[0];
    
    # turn method into sub
    s/^(\s*)method (.*)/$1sub $2 \#__METHOD/gm;
    
    return $_;
}

sub postfilter {
    $_ = $_[0];
    
    # turn sub back into method
    s/^(\s*)sub (.*?)\s* \#__METHOD/${1}method $2/gm;
    
    # the method->sub->method trick screws up comments that appear on the next
    # line; fix them now
    s/^(\s*method .+)\n\s+(\#.+?)\n(\s*)/$1\n$3$2\n$3/gm;
    
    # it also likes adding extra spaces between ; and # - remove these
    s/^(\s*\S[^#]+?)([;,]) +\#/$1$2 #/gm;
    
    # don't have completely empty lines between things; instead indent to the
    # correct level
    s/\n\n([\s]+)(\S.*)/\n$1\n$1$2/gm;
    
    return $_;
}

1;
