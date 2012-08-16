
=head1 NAME

VRPipe::Steps::chunk_genomic_region - a step

=head1 DESCRIPTION

Generate a chromosomal regions file, split according to chunk size parameter,
from either a fasta reference index file or a specific regions file

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::chunk_genomic_region with VRPipe::StepRole {
    method options_definition {
        return {
            genomic_region_file => VRPipe::StepOption->create(description => 'Genomic regions source file; absolute path to either fasta genomic reference index file (fai suffix), or specific regions file with columns CHR,FROM,TO'),
            chrom_list          => VRPipe::StepOption->create(description => 'Space separated list of chromosomes', optional => 1, default_value => '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y'),
            chunk_size          => VRPipe::StepOption->create(description => 'Number of base pairs to have in a chunk', optional => 1, default_value => 1_000_000),
            chunk_overlap       => VRPipe::StepOption->create(description => 'Chunk overlap size', optional => 1, default_value => 0),
            ploidy => VRPipe::StepOption->create(description => "File defining the ploidy to be used to call in different regions of the genome, eg {default=>2, X=>[{ from=>1, to=>60_000, M=>1 },{ from=>2_699_521, to=>154_931_043, M=>1 },],Y=>[{ from=>1, to=>59_373_566, M=>1, F=>0 }]}", optional => 1),
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $genomic_region_file = Path::Class::File->new($options->{genomic_region_file});
            $self->throw("genomic_region_file must be an absolute path") unless $genomic_region_file->is_absolute;
            
            my $chrom_list    = $options->{chrom_list};
            my $chunk_size    = $options->{chunk_size};
            my $chunk_overlap = $options->{chunk_overlap};
            my $ploidy        = $options->{ploidy};
            
            my $chunked_regions_file      = $self->output_file(output_key => 'chunked_regions_file', basename => $genomic_region_file->basename . qq[.$chunk_size.$chunk_overlap.txt], type => 'txt');
            my $chunked_regions_file_path = $chunked_regions_file->path;
            my $genomic_region_path       = $genomic_region_file->stringify;
            
            my $req = $self->new_requirements(memory => 50, time => 1);
            my $cmd = "use VRPipe::Steps::chunk_genomic_region; VRPipe::Steps::chunk_genomic_region->write_chunked_regions_file('$genomic_region_path', '$chunked_regions_file_path', '$chrom_list', '$chunk_size', '$chunk_overlap', '$ploidy');";
            
            $self->dispatch_vrpipecode($cmd, $req, { output_files => [$chunked_regions_file] });
        };
    }
    
    method outputs_definition {
        return { chunked_regions_file => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'chromosomal regions file split according to chunk size') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generate a chromosomal regions file, split according to chunk size, from fasta reference index file or specific regions file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method write_chunked_regions_file (ClassName|Object $self: Str|File $genomic_region_file, Str|File $chunked_regions_file_path, Str $chrom_list, PositiveInt $chunk_size, Int $chunk_overlap, Str|File $ploidy) {
        my $chroms = [split(' ', $chrom_list)];
        my $ploidy_regions;
        if ($ploidy) {
            $ploidy_regions = do $ploidy;
        }
        my $ploidy_default = $$ploidy_regions{default} || 2;
        
        my $reg_file = VRPipe::File->create(path => $genomic_region_file);
        my $rfh = $reg_file->openr;
        
        my @chr_regions;
        while (my $line = <$rfh>) {
            # genomic reference index file
            if ($genomic_region_file =~ /.fai$/) {
                my ($chr, $from, $to);
                for my $regex (@$chroms) {
                    if (!($line =~ /^($regex)\t(\d+)/i)) { next; }
                    $chr  = $1;
                    $from = 1;
                    $to   = $2;
                    last;
                }
                if (!defined $chr) { next; }
                
                if (!exists($$ploidy_regions{$chr})) {
                    push @chr_regions, { chr => $chr, from => $from, to => $to, female => $ploidy_default, male => $ploidy_default };
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
                        push @chr_regions, { chr => $chr, from => $from, to => $start - 1, female => $ploidy_default, male => $ploidy_default };
                    }
                    push @chr_regions, { chr => $chr, from => $start, to => $end, female => $female, male => $male };
                    $from = $end + 1;
                }
                
                if ($from < $to) {
                    push @chr_regions, { chr => $chr, from => $from, to => $to, female => $ploidy_default, male => $ploidy_default };
                }
            }
            # specific regions file
            else {
                chomp($line);
                if (!($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s*/)) { $self->throw("Could not parse $genomic_region_file: [$line]"); }
                push @chr_regions, { chr => $1, from => $2, to => $3, female => $ploidy_default, male => $ploidy_default };
            }
        }
        
        my $outfile = VRPipe::File->get(path => file($chunked_regions_file_path));
        my $ofh = $outfile->openw;
        
        my @chunks;
        for my $region (@chr_regions) {
            my $pos     = $$region{from};
            my $end_pos = $$region{to};
            my $female  = $$region{female};
            my $male    = $$region{male};
            while ($pos < $end_pos) {
                my $from = $pos;
                my $to   = $from + $chunk_size - 1;
                
                if ($to > $end_pos) { $to = $end_pos; }
                
                print $ofh join("\t", $$region{chr}, $from, $to, $female, "$male\n");
                
                $pos += $chunk_size - $chunk_overlap;
                if ($pos < 1) { $self->throw("The split size too small [$chunk_size]?\n"); }
            }
        }
    }
}

1;
