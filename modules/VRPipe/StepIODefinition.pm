=head1 NAME

VRPipe::StepIODefinition - define inputs and outputs of a Step

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

B<VRPipe> must know what a L<VRPipe::Step> requires as input, and what it
outputs. With this knowledge it can ensure that a B<VRPipe::Pipeline> makes
sense and that it is functioning correctly whilst running (it automatically
checks that input and output files exist before proceeding along the Pipeline).

A StepIODefinition can be used when writing a C<VRPipe::Steps::[step_name]>
module to define what is expected of input and output files without having to
know their paths (Steps should never care about file locations or names; they
are unknowable at the time of writing a Step anyway).

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

use VRPipe::Base;

class VRPipe::StepIODefinition extends VRPipe::Persistent {
    has 'type' => (is => 'rw',
                   isa => FileType,
                   coerce => 1,
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'min_files' => (is => 'rw',
                        isa => IntSQL[4],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1,
                        default => 1,
                        allow_key_to_default => 1);
    
    has 'max_files' => (is => 'rw',
                        isa => IntSQL[4],
                        traits => ['VRPipe::Persistent::Attributes'],
                        is_key => 1,
                        default => 1, # -1 means no maximum
                        allow_key_to_default => 1);
    
    has 'metadata' => (is => 'rw',
                       isa => 'HashRef',
                       traits => ['VRPipe::Persistent::Attributes'],
                       is_key => 1,
                       default => sub { {} },
                       allow_key_to_default => 1);
    
    has 'check_existence' => (is => 'rw',
                              isa => 'Bool',
                              traits => ['VRPipe::Persistent::Attributes'],
                              is_key => 1,
                              default => 1,
                              allow_key_to_default => 1);
    
    has 'description' => (is => 'rw',
                          isa => Varchar[255],
                          traits => ['VRPipe::Persistent::Attributes'],
                          is_key => 1);
    
    __PACKAGE__->make_persistent();
    
    method required_metadata_keys {
        my %metadata = %{$self->metadata}; # copy it so we don't alter original on next line, letting us reuse this same instance safely
        my %optional = map { $_ => 1 } @{delete $metadata{optional} || []};
        
        my @required;
        foreach my $key (sort keys %metadata) {
            next if exists $optional{$key};
            push(@required, $key);
        }
        
        return @required;
    }
}

1;