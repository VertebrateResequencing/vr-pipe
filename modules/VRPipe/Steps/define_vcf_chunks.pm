
=head1 NAME

VRPipe::Steps::define_vcf_chunks - a step

=head1 DESCRIPTION

This step defines chunk coordinates for input VCF files given a fixed number 
of sites per chunk and in buffer regions. The step outputs the chunk
coordinates  in bed format. Note: the step does not split the input VCF files
as it is  meant to precede steps which extract and process the VCF chunks on
the fly.

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk> and Petr Danecek <pd3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::define_vcf_chunks with VRPipe::StepRole {
    method options_definition {
        return {
            chunk_nsites => VRPipe::StepOption->create(
                description   => 'number of sites per chunk',
                optional      => 0,
                default_value => 100000
            ),
            buffer_nsites => VRPipe::StepOption->create(
                description   => 'number of sites in buffer regions',
                optional      => 1,
                default_value => 0
            ),
            max_chr_len => VRPipe::StepOption->create(
                description   => 'length of largest chromosome',
                optional      => 1,
                default_value => 249250621
            ),
            ref_vcf => VRPipe::StepOption->create(
                description => 'path to an indexed reference VCF which may contain the string "{CHROM}" that will be expanded to chromosome names',
                optional    => 1,
            ),
            chunk_by_ref => VRPipe::StepOption->create(
                description   => 'set to 1 to define chunks based on sites in the reference VCF (when multiple reference VCFs are given, only the first is used to define the chunks)',
                optional      => 1,
                default_value => 0
            ),
            regions => VRPipe::StepOption->create(
                description => 'regions to include (format should follow chr:start-end seperated by comma)',
                optional    => 1,
            ),
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools exe',
                optional      => 1,
                default_value => 'bcftools'
            ),
            tabix_exe => VRPipe::StepOption->create(
                description   => 'path to your tabix exe',
                optional      => 1,
                default_value => 'tabix'
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                min_files   => 1,
                max_files   => -1,
                description => '1 or more indexed VCF files',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $chunk_nsites  = $options->{chunk_nsites};
            my $buffer_nsites = $options->{buffer_nsites};
            my $ref_vcf       = $options->{ref_vcf};
            my $chunk_by_ref  = $options->{chunk_by_ref};
            my $max_chr_len   = $options->{max_chr_len};
            my $regions       = $options->{regions};
            my $bcftools_exe  = $options->{bcftools_exe};
            my $tabix_exe     = $options->{tabix_exe};
            
            my $vcf_to_chunk;
            if ($chunk_by_ref) {
                if ($ref_vcf) {
                    my @ref_vcfs = split(/\s+/, $ref_vcf);
                    $vcf_to_chunk = $ref_vcfs[0];
                    if ($vcf_to_chunk =~ /{CHROM}/ && !defined $regions) {
                        $self->throw("option regions is required when ref_vcf contains the string {CHROM}");
                    }
                }
                else {
                    $self->throw("chunk_by_ref is set to 1 while ref_vcf is not given!");
                }
            }
            
            foreach my $in_vcf (@{ $self->inputs->{vcf_files} }) {
                $vcf_to_chunk = $in_vcf->path->stringify unless (defined $vcf_to_chunk);
                my @regions = $self->define_regions($vcf_to_chunk, $regions, $tabix_exe);
                my @vcf_path = split(/\//, $vcf_to_chunk);
                pop @vcf_path;
                foreach my $region (@regions) {
                    my $chr = $region;
                    $chr =~ s/:.*$//;
                    my $outdir = join('/', @vcf_path) . "/$chr";
                    my $bed_file = $self->output_file(output_dir => "$outdir", output_key => 'bed_files', basename => "$region.${chunk_nsites}_$buffer_nsites.bed", type => 'txt');
                    my $bed_path = $bed_file->path;
                    my $this_cmd = "use VRPipe::Steps::define_vcf_chunks; VRPipe::Steps::define_vcf_chunks->define_chunks(outfile => q[$bed_path], in_path => q[$vcf_to_chunk], region => q[$region], buffer_nsites => q[$buffer_nsites], chunk_nsites => q[$chunk_nsites], max_chr_len => q[$max_chr_len], bcftools => q[$bcftools_exe]);";
                    my $req      = $self->new_requirements(memory => 1000, time => 1);
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$bed_file], block_and_skip_if_ok => 1 });
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            bed_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => -1,
                description => 'bed files containing coordinates of chunks',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Determine chunk coordinates for input VCF files - including a fixed number of sites per chunk";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method expand_chrom (ClassName|Object $self: Str $path, Str $region) {
        if (!defined $region) { return $path; }
        $region =~ s/:.*$//;
        $path =~ s/{CHROM}/$region/g;
        return $path;
    }
    
    method define_regions (Str $path, Str $regions, Str $tabix) {
        if ($regions) {
            my @regions = split(/\s*,\s*/, $regions);
            return @regions;
        }
        else {
            $self->throw("File $path must be tabix indexed!") unless -e "$path.tbi";
            my @chrs = grep { chomp } `$tabix -l $path`;
            return @chrs;
        }
    }
    
    method define_chunks (ClassName|Object $self: Str|File :$outfile!, Str :$in_path!, Str :$region!, Str :$buffer_nsites!, Str :$chunk_nsites!, Str :$max_chr_len!, Str :$bcftools!) {
        my $tot_sites = $buffer_nsites + $chunk_nsites;
        my (@chunks, @buffer);
        $in_path = $self->expand_chrom($in_path, $region);
        if ($region eq '.') { $region = ''; }
        my $cmd = "$bcftools view -g ^miss $in_path $region |";
        open(my $in, $cmd) or $self->throw("$cmd: $!");
        while (my $line = <$in>) {
            if (substr($line, 0, 1) eq '#') { next; }
            my $i = index($line, "\t");
            if ($i < 0) { $self->throw("Could not parse the line [CHR]: $line"); }
            my $chr = substr($line, 0, $i);
            my $j = index($line, "\t", $i + 1);
            if ($j < 0) { $self->throw("Could not parse the line [POS]: $line"); }
            my $pos = substr($line, $i + 1, $j - $i - 1);
            if (@buffer && $buffer[0][0] ne $chr or @buffer > $tot_sites) {
                my $chr_from = $buffer[0][0];
                my $pos_from = $buffer[0][1];
                my $pos_to   = $buffer[-1][1];
                my $nout     = @buffer;
                push @chunks, { chr => $chr_from, from => $pos_from, to => $pos_to, n => $nout };
                if ($chunk_nsites < @buffer) { splice(@buffer, 0, $chunk_nsites); }
                else                         { @buffer = (); }
            }
            push @buffer, [$chr, $pos];
        }
        if (@buffer) {
            my $chr_from = $buffer[0][0];
            my $pos_from = $buffer[0][1];
            my $pos_to   = $buffer[-1][1];
            my $nout     = @buffer;
            push @chunks, { chr => $chr_from, from => $pos_from, to => $pos_to, n => $nout };
        }
        close($in) or $self->throw("close $cmd");
        #if ( !@chunks ) { $self->throw("No chunks defined: $cmd\n"); }
        if (@chunks > 1 && $chunks[-1]{n} < $tot_sites * 0.75 && $chunks[-1]{chr} eq $chunks[-2]{chr}) {
            my $chunk = splice(@chunks, -1, 1);
            $chunks[-1]{to} = $$chunk{to};
            $chunks[-1]{n} += $$chunk{n};
        }
        
        if (!($region =~ /:/)) { # Whole genome or whole chromosome was requested. When on a new
            # chromosome, expand the first and last record to accompany
            # sites which may be present only in one (known_vcf vs in_vcf)
            for (my $i = 0; $i < @chunks; $i++) {
                if ($i == 0) { $chunks[0]{from} = 0; next; }
                if ($chunks[$i]{chr} ne $chunks[$i - 1]{chr}) {
                    $chunks[$i - 1]{to} = $max_chr_len; # last chunk, longest chr of human genome
                    $chunks[$i]{from} = 0;
                }
            }
            $chunks[-1]{to} = $max_chr_len;
        }
        
        my $out_file = VRPipe::File->get(path => $outfile);
        my $ofh = $out_file->openw;
        $out_file->disconnect;
        print $ofh "#chr\tstart\tend\tnsites\n";
        for my $chunk (@chunks) {
            print $ofh "$$chunk{chr}\t$$chunk{from}\t$$chunk{to}\t$$chunk{n}\n";
            # Here we could also run tabix to split the vcf file, however, for scalability,
            # each vcf chunk should ideally be generated using a separate VRPipe dispatch.
            # Therefore, vcf split, if required, should be carried out in a separate step.
        }
        $out_file->close;
        return 1;
    }

}

1;
