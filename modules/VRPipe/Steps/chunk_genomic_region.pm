use VRPipe::Base;

class VRPipe::Steps::chunk_genomic_region with VRPipe::StepRole {
	method options_definition {
		return {
				chunking_regions_file => VRPipe::StepOption->get( description => 'target regions source; absolute path to either fasta genomic ref index file (fai suffix), or specific regions file with columns CHR,FROM,TO'),
				chrom_list =>  VRPipe::StepOption->get(description => 'Names of Chromosomes to split', optional => 1, default_value => '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y' ),
				chunk_size => VRPipe::StepOption->get(description => 'Split file chunk size', optional => 1, default_value => 100000000),
				chunks_overlap => VRPipe::StepOption->get(description => 'chunk overlap size', optional => 1, default_value => 0),
				pseudo_autosomal_definition => VRPipe::StepOption->get(description => "structure definition if restricting sex chromosomes to pseudo_autosomal regions, eg {X=>[{region=>'1-60000'},{region=>'2699521-154931043'}],Y=>[{region=>'1-59373566'},],}", optional => 1),
		};
	}
    method inputs_definition {
        return { };
    }

    method body_sub {
		return sub {

			my $self = shift;
			my $options = $self->options;
			my $chunking_regions_file = $options->{chunking_regions_file};
			my $chrom_list = $options->{chrom_list};
			my $chunk_size = $options->{chunk_size};
			my $chunks_overlap = $options->{chunks_overlap};
			my $pseudo_autosomal_definition = $options->{pseudo_autosomal_definition};
			$pseudo_autosomal_definition  =~ s/'/#/g;	# avoids parsing problems when passing args

			my $chunked_regions_file = $self->output_file(output_key => 'chunked_regions_file', basename => 'regions.txt', type => 'txt');
			my $chunked_regions_file_path = $chunked_regions_file->path;

            my $req = $self->new_requirements(memory => 500, time => 1);
			my $cmd = "use VRPipe::Steps::chunk_genomic_region; VRPipe::Steps::chunk_genomic_region->write_chunked_regions_file('$chunking_regions_file', '$chunked_regions_file_path', '$chrom_list', '$chunk_size', '$chunks_overlap', '$pseudo_autosomal_definition');";

			$self->dispatch_vrpipecode($cmd, $req, {output_files => [$chunked_regions_file]});
		};
    }

    method outputs_definition {
        return { chunked_regions_file => VRPipe::StepIODefinition->get(type => 'txt', max_files => 1, description => 'chromosomal regions file split according to chunk size') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generate a chromosomal regions file, split according to chunk size, from fasta reference file or specific regions file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }

    method write_chunked_regions_file (ClassName|Object $self: Str|File $chunking_regions_file, Str|File $chunked_regions_file_path, Str $chrom_list, PositiveInt $chunk_size, Int $chunks_overlap, Str $pseudo_autosomal_definition) {
        
			my $chroms = [split(' ',$chrom_list)];
			my $pseudo_autosomal;
			if ($pseudo_autosomal_definition) {
				$pseudo_autosomal_definition  =~ s/#/'/g;
				eval '$pseudo_autosomal='.$pseudo_autosomal_definition;
				die "Failed to define pseudo_autosomal:$@" if $@;
			}

			my $reg_file = VRPipe::File->get(path => $chunking_regions_file);
			my $rfh = $reg_file->openr;

			my @chr_regions;
			while (my $line=<$rfh>) {
				if ($chunking_regions_file =~ /.fai$/) { # genomic reference index file
					my ($chr,$from,$to);
					for my $regex (@$chroms)
					{
						if ( !($line=~/^($regex)\t(\d+)/i) ) { next; }
						$chr  = $1;
						$from = 1;
						$to   = $2;
						last;
					}
					if ( !defined $chr ) { next; }

					if ( !exists($pseudo_autosomal->{$chr}) )
					{
						push @chr_regions, { chr=>$chr, from=>$from, to=>$to };
						next;
					}
					for my $reg (@{$pseudo_autosomal->{$chr}})
					{
						my ($start,$end) = split(/-/,$$reg{region});
						if ( $start>$from )
						{
							push @chr_regions, { chr=>$chr, from=>$from, to=>$start-1 };
						}
						push @chr_regions, { chr=>$chr, from=>$start, to=>$end };
						$from = $end+1;
					}

					if ( $from<$to ) {
						push @chr_regions, { chr=>$chr, from=>$from, to=>$to };
					}
				}
				else {	# specific regions file
					chomp($line);
					if ( !($line=~/^(\S+)\s+(\d+)\s+(\d+)\s*$/) ) { $self->throw("Could not parse $chunking_regions_file: [$line]"); }
					push @chr_regions, { chr=>$1, from=>$2, to=>$3 };
				}
			}

            my $outfile = VRPipe::File->get(path => file($chunked_regions_file_path));
			my $ofh = $outfile->openw;

			my @chunks;
			for my $region (@chr_regions)
			{
				my $pos     = $$region{from};
				my $end_pos = $$region{to};
				while ($pos<$end_pos) {
					my $from = $pos;
					my $to   = $from+$chunk_size-1;
	
					if ( $to>$end_pos ) { $to=$end_pos; }
	
					print $ofh join("\t", $$region{chr}, $from, "$to\n");
	
					$pos += $chunk_size - $chunks_overlap;
					if ( $pos<1 ) { $self->throw("The split size too small [$chunk_size]?\n"); }
				}
			}

    } # write_chunked_regions_file

}

1;
