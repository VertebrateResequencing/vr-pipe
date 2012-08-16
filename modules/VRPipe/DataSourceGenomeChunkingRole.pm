
=head1 NAME

VRPipe::DataSourceGenomeChunkingRole - role extends datasources to be  split
into chunks across the genome

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

role VRPipe::DataSourceGenomeChunkingRole with VRPipe::DataSourceRole {
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
    
    has 'chunk_override_file' => (
        is      => 'rw',
        isa     => 'Str',
        lazy    => 1,
        builder => 'options'
    );
    
    has 'chunk_size' => (
        is      => 'rw',
        isa     => 'Int',
        lazy    => 1,
        default => 1000000
    );
    
    has 'chunk_overlap' => (
        is      => 'rw',
        isa     => 'Int',
        lazy    => 1,
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
    
    around method_description (Str $method) {
        my $description = $self->$orig($method);
        return $description . q[ Each dataelement will be duplicated in chunks across the genome. The option 'reference_index' is the absolute path to the fasta index (.fai) file associated with the reference fasta file, 'chunk_override_file' is a file defining chunk specific options that may be overridden (required, but may point to an empty file), 'chunk_size' the size of the chunks in bp, 'chunk_overlap' defines how much overlap to have beteen chunks, 'chrom_list' (a space separated list) will restrict to specified the chromosomes (must match chromosome names in dict file), 'ploidy' is an optional file specifying the ploidy to be used for males and females in defined regions of the genome, eg {default=>2, X=>[{ from=>1, to=>60_000, M=>1 },{ from=>2_699_521, to=>154_931_043, M=>1 },],Y=>[{ from=>1, to=>59_373_566, M=>1, F=>0 }]}.];
    }
    
    around _create_elements (ArrayRef $e_args!) {
        my $chunks = $self->chunks();
        my @element_args;
        foreach my $e_arg (@$e_args) {
            my $result       = $e_arg->{result};
            my $meta         = $self->common_metadata($result->{paths});
            my $chrom_select = $meta->{chrom} || $meta->{split_sequence} || '';
            $chrom_select =~ s/^chrom//;
            foreach my $chunk (@$chunks) {
                if ($chrom_select) {
                    my $chrom = $chunk->{chrom};
                    $chrom =~ s/^chrom//;
                    next if ($chrom ne $chrom_select);
                }
                $result = { %$result, %$chunk, chunk_override_file => $self->chunk_override_file };
                push @element_args, { %$e_arg, result => $result };
            }
        }
        $self->$orig(\@element_args);
    }
    
    method common_metadata (ArrayRef[Str] $file_paths!) {
        my %meta;
        foreach my $file_path (@$file_paths) {
            my $file = VRPipe::File->get(path => $file_path);
            my $file_meta = $file->metadata;
            foreach my $key (keys %$file_meta) {
                $meta{$key}->{ $$file_meta{$key} } += 1;
            }
        }
        # Only keep metadata common to all files
        my $common_meta = {};
        foreach my $key (keys %meta) {
            my @vals = keys %{ $meta{$key} };
            next unless (@vals == 1 && $meta{$key}->{ $vals[0] } == @$file_paths);
            $common_meta->{$key} = $vals[0];
        }
        return $common_meta;
    }
    
    around options {
        my $options = $self->$orig();
        $self->reference_index(delete $options->{reference_index});
        $self->chunk_override_file(delete $options->{chunk_override_file});
        $self->chunk_size(delete $options->{chunk_size})       if (exists $options->{chunk_size});
        $self->chunk_overlap(delete $options->{chunk_overlap}) if (exists $options->{chunk_overlap});
        if (exists $options->{chrom_list}) {
            my $chrom_list = delete $options->{chrom_list};
            $self->chrom_list([split(' ', $chrom_list)]) if $chrom_list;
        }
        if (exists $options->{ploidy}) {
            my $ploidy = delete $options->{ploidy};
            $self->ploidy(do $ploidy) if $ploidy;
        }
        return $options;
    }
    
    around method_options (Str $method) {
        my @return = ();
        push @return, $self->$orig("$method");
        push @return, $self->$orig('_with_genome_chunking');
        return @return;
    }
    
    method _with_genome_chunking (Defined :$handle!, Str|File :$reference_index!, Str|File :$chunk_override_file!, Int :$chunk_size! = 1000000, Int :$chunk_overlap! = 0, Str :$chrom_list?, Str|File :$ploidy?) {
        return;
    }
    
    method _build_chunks {
        my $reference_index = $self->reference_index;
        my $chunk_size      = $self->chunk_size;
        my $chunk_overlap   = $self->chunk_overlap;
        my $chrom_list      = $self->chrom_list;
        my $ploidy_regions  = $self->ploidy;
        
        my $ploidy_default = $$ploidy_regions{default} || 2;
        
        my $pars = VRPipe::Parser->create('fai', { file => $reference_index });
        my $parsed_record = $pars->parsed_record();
        
        my @chroms = @$chrom_list;
        unless (@chroms) {
            @chroms = $pars->seq_ids;
        }
        
        my @chr_regions;
        while ($pars->next_record()) {
            my ($chr, $from, $to);
            for my $regex (@chroms) {
                my $seq = $parsed_record->{SN};
                if (!($seq =~ /^([cC]hr|[cC]hrom)?($regex)$/i)) { next; }
                $chr  = $seq;
                $from = 1;
                $to   = $parsed_record->{LN};
                last;
            }
            if (!defined $chr) { next; }
            
            if (!exists($$ploidy_regions{$chr})) {
                push @chr_regions, { chrom => $chr, from => $from, to => $to, female => $ploidy_default, male => $ploidy_default };
                next;
            }
            
            # Split the chunks as necessary
            if (exists $$ploidy_regions{$chr} && ref($$ploidy_regions{$chr}) ne 'ARRAY') { $self->throw("Ploidy definition not well defined for $chr regions.\n"); }
            foreach my $reg (sort { $$a{from} <=> $$b{to} } @{ $$ploidy_regions{$chr} }) {
                my $start  = $$reg{from};
                my $end    = $$reg{to};
                my $female = defined $$reg{F} ? $$reg{F} : $ploidy_default;
                my $male   = defined $$reg{M} ? $$reg{M} : $ploidy_default;
                if ($start > $from) {
                    push @chr_regions, { chrom => $chr, from => $from, to => $start - 1, female => $ploidy_default, male => $ploidy_default };
                }
                push @chr_regions, { chrom => $chr, from => $start, to => $end, female => $female, male => $male };
                $from = $end + 1;
            }
            
            if ($from < $to) {
                push @chr_regions, { chrom => $chr, from => $from, to => $to, female => $ploidy_default, male => $ploidy_default };
            }
        }
        
        my @chunks;
        my $seq_no = 1;
        for my $region (@chr_regions) {
            my $pos     = $$region{from};
            my $end_pos = $$region{to};
            my $female  = $$region{female};
            my $male    = $$region{male};
            while ($pos < $end_pos) {
                my $from = $pos;
                my $to   = $from + $chunk_size - 1;
                
                if ($to > $end_pos) { $to = $end_pos; }
                
                if (keys %$ploidy_regions) {
                    push @chunks, { chrom => $$region{chrom}, from => $from, to => $to, female_ploidy => $female, male_ploidy => $male, seq_no => $seq_no };
                }
                else {
                    push @chunks, { chrom => $$region{chrom}, from => $from, to => $to, seq_no => $seq_no };
                }
                ++$seq_no;
                
                $pos += $chunk_size - $chunk_overlap;
                if ($pos < 1) { $self->throw("The split size too small [$chunk_size]?\n"); }
            }
        }
        return \@chunks;
    }
}

1;
