=head1 NAME

VRPipe::Base::ReadWritable - a role for objects that represent a file

=head1 SYNOPSIS

use MooseX::Declare;
class VRPipe::MyClass with (VRPipe::Base::ReadWritable) {
    #...
    
    method read () {
        # ...
    }
    method write () {
        # ...
    }
}

=head1 DESCRIPTION

This role describes the basics of a file: it is something that string data can
be read from, and that string data can be written to.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

role VRPipe::Base::ReadWritable {
    use MooseX::Aliases;
    
=head2 read

 Title   : read (alias next_line)
 Usage   : my $line = $obj->read();
 Function: Read the next line from the file.
 Returns : string
 Args    : n/a

=cut
    requires 'read';
    alias next_line => 'read';
    
=head2 write

 Title   : write (alias print)
 Usage   : $obj->write("string");
 Function: Output the given string to the file.
 Returns : n/a
 Args    : string

=cut  
    requires 'write';
    alias print => 'write';
}
