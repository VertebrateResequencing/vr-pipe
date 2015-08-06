
=head1 NAME

VRPipe::StepGenomeChunkingRole - role extends steps to be split into chunks
across the genome

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

role VRPipe::StepGenomeChunkingRole with VRPipe::StepRole {
    has 'chunks' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        lazy    => 1,
        builder => '_build_chunks'
    );
    
    has 'reference_index' => (
        is      => 'rw',
        isa     => 'Str',
        lazy    => 1,
        builder => 'options'
    );
    
    has 'chunk_size' => (
        is      => 'rw',
        isa     => 'Int',
        default => 50000000
    );
    
    has 'chunk_overlap' => (
        is      => 'rw',
        isa     => 'Int',
        default => 0
    );
    
    has 'chrom_list' => (
        is      => 'rw',
        isa     => 'ArrayRef',
        lazy    => 1,
        default => sub { [] }
    );
    
    has 'ploidy' => (
        is      => 'rw',
        isa     => 'HashRef',
        lazy    => 1,
        default => sub { {} }
    );
    
    has 'target_regions' => (
        is      => 'rw',
        isa     => 'Str',
        lazy    => 1,
        default => ''
    );
    
    has '_processed_extra_options' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    around options {
        my $options = $self->$orig();
        unless ($self->_processed_extra_options) {
            $self->reference_index(file($options->{reference_fasta} . '.fai')->absolute->stringify);
            $self->chunk_size(delete $options->{chunk_size})         if (exists $options->{chunk_size});
            $self->chunk_overlap(delete $options->{chunk_overlap})   if (exists $options->{chunk_overlap});
            $self->target_regions(delete $options->{target_regions}) if (exists $options->{target_regions});
            if (exists $options->{chrom_list}) {
                my $chrom_list = delete $options->{chrom_list};
                $self->chrom_list([split(' ', $chrom_list)]) if $chrom_list;
            }
            if (exists $options->{ploidy}) {
                my $ploidy = delete $options->{ploidy};
                $self->ploidy(do $ploidy) if $ploidy;
            }
            $self->_processed_extra_options(1);
        }
        return $options;
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            reference_fasta => VRPipe::StepOption->create(description => 'Absolute path to the reference fasta file. The fasta index (fai) file must exist'),
            chrom_list      => VRPipe::StepOption->create(description => 'Space separated list of chromosomes', optional => 1),
            chunk_size      => VRPipe::StepOption->create(description => 'Number of base pairs to have in a chunk', optional => 1, default_value => 50_000_000),
            chunk_overlap   => VRPipe::StepOption->create(description => 'Chunk overlap size', optional => 1, default_value => 0),
            target_regions => VRPipe::StepOption->create(description => "BED file defining targetted sequecing regions. Chunks will be created such that there is *at least* chunk_size bases in a chunk made from ",                                                                                   optional => 1),
            ploidy         => VRPipe::StepOption->create(description => "File defining the ploidy to be used to call in different regions of the genome, eg {default=>2, X=>[{ from=>1, to=>60_000, M=>1 },{ from=>2_699_521, to=>154_931_043, M=>1 },],Y=>[{ from=>1, to=>59_373_566, M=>1, F=>0 }]}", optional => 1),
        };
    }
    
    method _build_chunks {
        my $reference_index = $self->reference_index;
        my $chunk_size      = $self->chunk_size;
        my $chunk_overlap   = $self->chunk_overlap;
        my $chrom_list      = $self->chrom_list;
        my $ploidy_regions  = $self->ploidy;
        my $target_regions  = $self->target_regions;
        
        use VRPipe::Utils::GenomeChunking;
        my $chunk_util = VRPipe::Utils::GenomeChunking->new();
        return $chunk_util->chunks(reference_index => $reference_index, chunk_size => $chunk_size, chroms => $chrom_list, ploidy => $ploidy_regions, $target_regions ? (target_regions => $target_regions) : ());
    }
}

1;
