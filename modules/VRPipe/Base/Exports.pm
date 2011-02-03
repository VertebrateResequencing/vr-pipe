=head1 NAME

VRPipe::Base::Exports - commonly used methods from external CPAN modules

=head1 SYNOPSIS

package Foo;
use Moose;


=head1 DESCRIPTION

Provides a grab-bag of very commonly used methods from external CPAN modules, so
you can call $self->common_method instead of { use CPAN::module;
CPAN::module->common_method }.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

role VRPipe::Base::Exports {
    use MooseX::StrictConstructor;
    use Cwd qw(abs_path cwd);
    use File::Spec;
    use File::Basename;
    require File::Path;
    require File::Copy;
    
=head2 catfile

 Title   : catfile
 Usage   : my ($path) = $obj->catfile('dir', 'subdir', 'filename'); 
 Function: Constructs a full pathname in a cross-platform safe way. Just an
           alias to File::Spec->catfile.
 Returns : the full path
 Args    : as per File::Spec->catfile

=cut
    method catfile (@parts) {
        return File::Spec->catfile(@parts);
    }
}

1;
