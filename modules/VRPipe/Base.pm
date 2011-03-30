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

        # if you do something that will require special consideration when
        # destroying ourselves, write a method to handle that and register it
        $self->register_for_cleanup('my_cleanup_method');
    }
}

=head1 DESCRIPTION

Used in place of 'use MooseX::Declare;', providing the same functionality. In
addition it sets the base class of your class to VRPipe::Base:UseMoose, which in
turn gives you the VRPipe::Base::Debuggable and Cleanable roles.

Furthermore, all VRPipe::Base::Types are loaded for you, and
MooseX::StrictConstructor is turned on.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Base extends MooseX::Declare is dirty {
    use aliased 'VRPipe::Base::Declare::Syntax::Keyword::Class', 'ClassKeyword';
    use aliased 'VRPipe::Base::Declare::Syntax::Keyword::Role',  'RoleKeyword';
    use aliased 'MooseX::Declare::Syntax::Keyword::Namespace',   'NamespaceKeyword';
    
    clean;
    
    around keywords (ClassName $self:) {
        ClassKeyword->new(identifier => 'class'),
        RoleKeyword->new(identifier => 'role'),
        NamespaceKeyword->new(identifier => 'namespace'),
    }
}

1;