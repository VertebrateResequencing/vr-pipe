=head1 NAME

VRPipe::Base::FHWrapper - treat a file handle like an object

=head1 SYNOPSIS

my $obj = VRPipe::Base::FHWrapper->new(handle => $fh);
my $string = $obj->read;
$obj->write($string);

=head1 DESCRIPTION

Used when coercing a file handle to the VRPipe::Base::Types 'File' type.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::Base::FHWrapper with (VRPipe::Base::ReadWritable) {
    has 'handle' => ( is  => 'rw',
                      isa => 'FileHandle',
                      required => 1 );
    
    method read {
        my $fh = $self->handle();
        return <$fh>;
    }
    
    method write {
        my $fh = $self->handle();
        print {$fh} @_;
    }
}
