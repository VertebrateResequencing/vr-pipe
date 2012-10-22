
=head1 NAME

VRPipe::Steps::vcf_split - a step

=head1 DESCRIPTION

Uses tabix to split indexed, compressed VCF files into genomic regions

=head1 AUTHOR

Chris Joyce    <cj5@sanger.ac.uk>. Shane McCarthy <sm15@sanger.ac.uk>.

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
        return {
            tabix_exe       => VRPipe::StepOption->create(description => 'path to tabix executable', optional => 1, default_value => 'tabix'),
            reference_index => VRPipe::StepOption->create(description => 'absolute path to the fasta index (.fai) file associated with reference fasta file'),
            chunk_size => VRPipe::StepOption->create(description => 'Chunk size to split the vcf into.', default_value => 10000000)
        };
    }
    
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', description => '1 or more indexed, compressed vcf files to split', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self       = shift;
            my $options    = $self->options;
            my $tabix      = $options->{tabix_exe};
            my $chunk_size = $options->{chunk_size};
            
            my $fai_file = file($options->{reference_index});
            $self->throw("reference_index must be an absolute path") unless $fai_file->is_absolute;
            
            my $pars = VRPipe::Parser->create('fai', { file => $fai_file });
            my %seq_lengths = $pars->seq_lengths;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $vcf_meta = $vcf->metadata;
                my $seq_no   = 0;
                my $vcf_path = $vcf->path;
                open(my $th, "$tabix -l $vcf_path |") or $self->throw("$tabix -l $vcf_path: $!");
                while (my $chr = <$th>) {
                    chomp($chr);
                    my $from = 1;
                    while ($from < $seq_lengths{$chr}) {
                        my $to = $from + $chunk_size - 1;
                        if ($to > $seq_lengths{$chr}) { $to = $seq_lengths{$chr}; }
                        
                        my $basename       = "${chr}_${from}-${to}." . $vcf->basename;
                        my $vcf_split_file = $self->output_file(
                            output_key => 'vcf_split_files',
                            basename   => $basename,
                            type       => 'vcf',
                            metadata   => { %$vcf_meta, chrom => $chr, from => $from, to => $to, seq_no => $seq_no }
                        );
                        $seq_no++;
                        
                        my $vcf_split_path = $vcf_split_file->path;
                        
                        my $cmd = qq[$tabix -h $vcf_path $chr:$from-$to | bgzip -c > $vcf_split_path];
                        $self->dispatch([$cmd, $req, { output_files => [$vcf_split_file] }]);
                        $from = $to + 1;
                    }
                }
                close($th);
            }
        };
    }
    
    method outputs_definition {
        return {
            vcf_split_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'a set of split vcf.gz files for each input vcf',
                metadata    => { seq_no => 'a sequence number assigned by the split for reassembly in correct order' }
            )
        };
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
