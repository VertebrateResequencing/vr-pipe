
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
    
    has 'result' => (
        is     => 'rw',
        isa    => 'HashRef',
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'withdrawn' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    __PACKAGE__->make_persistent(); # has_many => [element_states => 'VRPipe::DataElementState'] doesn't work because of ordering issues?
    
    method element_states {
        return VRPipe::DataElementState->search({ dataelement => $self->id }, { prefetch => [qw(dataelement pipelinesetup)] });
    }
    
    method files {
        my $paid = $self->result->{paths} || return;
        if (ref($paid)) {
            # we used to store an array ref of file path strings in the db
            return [map { VRPipe::File->get(path => $_) } @$paid];
        }
        return [VRPipe::PersistentArray->get(id => $paid)->member_instances];
    }
    
    method paths {
        return map { $_->path->stringify } @{ $self->files || return };
    }
    
    method start_from_scratch (VRPipe::PipelineSetup $setup, ArrayRef[PositiveInt] $step_numbers?) {
        VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $self)->start_from_scratch($step_numbers ? $step_numbers : ());
    }
    
    # most datasources have {result}->{paths} which contains absolute file path
    # strings. We convert these to VRPipe::File objects, and then to a
    # PersistentArray and store just the id. This saves db space and ensures we
    # can match the result regardless of order of paths, reducing spurious
    # DataElement withdrawal and recreation
    around bulk_create_or_update (ClassName|Object $self: @args) {
        foreach my $args (@args) {
            $self->_deflate_paths($args->{result});
        }
        
        return $self->$orig(@args);
    }
    
    around _get (ClassName|Object $self: Bool $create, %args) {
        if (defined $args{result}) {
            $self->_deflate_paths($args{result});
        }
        
        return $self->$orig($create, %args);
    }
    
    around search_rs (ClassName|Object $self: HashRef $search_args, Maybe[HashRef] $search_attributes?) {
        if (defined $search_args->{result}) {
            $self->_deflate_paths($search_args->{result});
        }
        
        my @args = ($search_args);
        push(@args, $search_attributes) if $search_attributes;
        return $self->$orig(@args);
    }
    
    sub _deflate_paths {
        my ($self, $result) = @_;
        if (defined $result->{paths} && ref($result->{paths})) {
            $result->{paths} = VRPipe::PersistentArray->get(members => [map { VRPipe::File->get(path => $_) } @{ $result->{paths} }], any_order => 1)->id;
        }
    }
}

1;
