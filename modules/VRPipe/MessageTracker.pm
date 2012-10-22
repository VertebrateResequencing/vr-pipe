
=head1 NAME

VRPipe::MessageTracker - track messages that have been sent to avoid repetition

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

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

use VRPipe::Base;

class VRPipe::MessageTracker extends VRPipe::Persistent {
    use DateTime;
    use POSIX qw(ceil);
    
    has 'subject' => (
        is     => 'rw',
        isa    => Varchar [32],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'message' => (
        is      => 'rw',
        isa     => Varchar [32],
        traits  => ['VRPipe::Persistent::Attributes'],
        default => ''
    );
    
    __PACKAGE__->make_persistent();
    
    around get (ClassName|Object $self: %args) {
        if (defined $args{subject} && length($args{subject}) != 32) {
            $args{subject} = $self->digest($args{subject});
        }
        return $self->$orig(%args);
    }
    
    around create (ClassName|Object $self: %args) {
        if (defined $args{subject} && length($args{subject}) != 32) {
            $args{subject} = $self->digest($args{subject});
        }
        return $self->$orig(%args);
    }
    
    method digest (ClassName|Object $self: Str $str) {
        my $dmd5 = Digest::MD5->new();
        $dmd5->add($str);
        return $dmd5->hexdigest;
    }
    
    method already_sent (Str $message) {
        my $digest  = $self->digest($message);
        my $current = $self->message;
        if ($digest eq $current) {
            return 1;
        }
        else {
            $self->message($digest);
            $self->update;
            return 0;
        }
    }
}

1;
