
=head1 NAME

VRPipe::DataElement - an element of a DataSource

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A DataElement respresents an element of a DataSource, and is the input to one
or more pipeline steps.

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

class VRPipe::DataElement extends VRPipe::Persistent {
    has 'datasource' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::DataSource'
    );
    
    has 'filelist' => (
        is         => 'rw',
        isa        => Persistent,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::FileList'
    );
    
    has 'keyvallist' => (
        is         => 'rw',
        isa        => Persistent,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::KeyValList'
    );
    
    has 'withdrawn' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    our $last_paths;
    our $last_filelist;
    
    __PACKAGE__->make_persistent(); # has_many => [element_states => 'VRPipe::DataElementState'] doesn't work because of ordering issues?
    
    method element_states {
        return VRPipe::DataElementState->search({ dataelement => $self->id }, { prefetch => [qw(dataelement pipelinesetup)] });
    }
    
    method files {
        my $filelist = $self->filelist;
        return [$filelist->files];
    }
    
    method paths {
        return map { $_->path->stringify } @{ $self->files || return };
    }
    
    method metadata {
        my $keyvallist = $self->keyvallist;
        return $keyvallist->as_hashref;
    }
    
    method result {
        my $result = {};
        my @paths  = $self->paths;
        if (@paths) {
            $result->{paths} = \@paths;
        }
        my $meta = $self->metadata;
        while (my ($key, $val) = each %{$meta}) {
            $result->{$key} = $val;
        }
        return $result;
    }
    
    method start_from_scratch (VRPipe::PipelineSetup $setup, ArrayRef[PositiveInt] $step_numbers?) {
        VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $self)->start_from_scratch($step_numbers ? (step_numbers => $step_numbers) : ());
    }
    
    # most datasources have {result}->{paths} which contains absolute file path
    # strings. We convert these to VRPipe::File objects, and then to a
    # FileList and store just the id. This saves db space and ensures we
    # can match the result regardless of order of paths, preventing spurious
    # DataElement withdrawal and recreation
    around bulk_create_or_update (ClassName|Object $self: @args) {
        foreach my $args (@args) {
            $self->_convert_result($args);
        }
        
        return $self->$orig(@args);
    }
    
    around _get (ClassName|Object $self: Bool $create, %args) {
        $self->_convert_result(\%args);
        return $self->$orig($create, %args);
    }
    
    around search_rs (ClassName|Object $self: HashRef $search_args, Maybe[HashRef] $search_attributes?) {
        $self->_convert_result($search_args);
        
        my @args = ($search_args);
        push(@args, $search_attributes) if $search_attributes;
        return $self->$orig(@args);
    }
    
    sub _convert_result {
        my ($self, $args) = @_;
        return unless exists $args->{result};
        my $result = delete $args->{result};
        
        # convert paths in the result to a filelist
        my $paths    = $result->{paths};
        my $protocol = delete $result->{protocol};
        if (defined $paths && ref($paths)) {
            # we might see the same set of $paths in multiple sequential calls
            # to _convert_result, and if there are thousands of paths we can
            # potentially save hours by caching the result
            my $paths_str = join('.', @$paths);
            if ($last_paths && $last_paths eq $paths_str) {
                $args->{filelist} = $last_filelist;
            }
            else {
                $args->{filelist} = VRPipe::FileList->get(files => [map { VRPipe::File->get(path => $_, $protocol ? (protocol => $protocol) : ()) } @{$paths}])->id;
                $last_paths       = $paths_str;
                $last_filelist    = $args->{filelist};
            }
        }
        $args->{filelist} ||= VRPipe::FileList->get(files => [])->id;
        
        # convert anything remaining in the result to a keyvals
        my $hash = {};
        while (my ($key, $val) = each %{$result}) {
            next if $key eq 'paths';
            $hash->{$key} = $val;
        }
        $args->{keyvallist} = VRPipe::KeyValList->get(hash => $hash)->id;
    }
}

1;
