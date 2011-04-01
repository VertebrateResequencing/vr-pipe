=head1 NAME

VRPipe::Base::FilePathWrapper - treat a file path like an object

=head1 SYNOPSIS

my $obj = VRPipe::Base::FilePathWrapper->new(path => "/path/to/file");
my $string = $obj->read;
$obj->write($string);

=head1 DESCRIPTION

Used when coercing a file path string to the VRPipe::Base::Types 'File' type.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Base::FilePathWrapper with (VRPipe::Base::ReadWritable) {
    has 'path' =>   ( is  => 'rw',
                      isa => 'Str',
                      required => 1 );
    has 'handle' => ( is => 'ro',
                      isa => 'FileHandle',
                      writer => '_set_handle' );
    
    method read {
        my $fh = $self->handle();
        return <$fh>;
    }
    
    method write {
        my $fh = $self->handle();
        print {$fh} @_;
    }
}
