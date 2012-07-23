
=head1 NAME

VRPipe::Steps::vcf_split - a step

=head1 DESCRIPTION

Uses tabix to split indexed, compressed VCF files into genomic regions, as defied by a file of chomosome region chunks
[we could implement other split strategies]

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::vcf_split with VRPipe::StepRole {
    method options_definition {
        return { tabix_exe => VRPipe::StepOption->create(description => 'path to tabix executable', optional => 1, default_value => 'tabix'), };
    }
    
    method inputs_definition {
        return { vcf_files            => VRPipe::StepIODefinition->create(type => 'vcf', description => '1 or more indexed, compressed vcf files to split', max_files   => -1),
                 chunked_regions_file => VRPipe::StepIODefinition->create(type => 'txt', max_files   => 1,                                                  description => 'file of chomosome region chunks into which each vcf should be split'), };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $tabix   = $options->{tabix_exe};
            my $req     = $self->new_requirements(memory => 500, time => 1);
            
            # Split based upon chunked_regions_file; we could implement other split strategies
            my %chunks;      # stores the chunk regions by chromosome
            my $chunk_file = $self->inputs->{chunked_regions_file}[0];
            my $cfh        = $chunk_file->openr;
            while (<$cfh>) {
                my ($chr, $from, $to) = split;
                push(@{ $chunks{$chr} }, "${from}-${to}");
            }
            
            my $seq_no = 0;
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                # Get chromosomes in vcf
                my $vcf_path = $vcf->path;
                open(CMD, "$tabix -l $vcf_path|") or $self->throw("$tabix -l $vcf_path: $!");
                while (my $chr = <CMD>) {
                    chomp($chr);
                    foreach my $region (@{ $chunks{$chr} }) {
                        my $basename = $vcf->basename;
                        $basename =~ s/\.vcf.gz$/.${chr}-${region}.vcf.gz/;
                        my $vcf_split_file = $self->output_file(output_key => 'vcf_split_files', basename => $basename, type => 'vcf', metadata => { seq_no => $seq_no });
                        $seq_no++;
                        
                        my $vcf_split_path = $vcf_split_file->path;
                        
                        my $cmd = qq[$tabix -h $vcf_path $chr:$region | bgzip -c > $vcf_split_path];
                        $self->dispatch([$cmd, $req, { output_files => [$vcf_split_file] }]);
                    }
                }
            }
        };
    }
    
    method outputs_definition {
        return { vcf_split_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a set of split vcf.gz files for each input vcf') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Split indexed, compressed vcf files by genomic regions";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
