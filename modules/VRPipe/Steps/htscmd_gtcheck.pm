
=head1 NAME

VRPipe::Steps::htscmd_gtcheck - a step

=head1 DESCRIPTION

Runs htscmd gtcheck to check sample identities in query vcf files against
genotypes in vcf file.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::htscmd_gtcheck with VRPipe::StepRole {
    method options_definition {
        return {
            expected_sample_from_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key that expected_sample should be taken from',
                optional      => 1,
                default_value => 'sample'
            ),
            genotypes_vcf => VRPipe::StepOption->create(description => 'absolute path to genotypes vcf file'),
            htscmd_exe    => VRPipe::StepOption->create(
                description   => 'path to the htscmd executable',
                optional      => 1,
                default_value => 'htscmd'
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'vcf files for genotyping',
                metadata    => {
                    source_bam => 'input bam path',
                    sample     => 'sample name'
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $genotypes_vcf = file($options->{genotypes_vcf});
            $self->throw("genotypes_vcf must be an absolute path") unless $genotypes_vcf->is_absolute;
            my $htscmd_exe   = $options->{htscmd_exe};
            my $expected_key = $options->{expected_sample_from_metadata_key};
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $vcf_path   = $vcf->path;
                my $meta       = $vcf->metadata;
                my $sample     = $meta->{$expected_key};
                my $source_bam = $meta->{source_bam};
                
                my $gtypex_file = $self->output_file(
                    output_key => 'htscmd_gtcheck_files',
                    basename   => $vcf->basename . '.gtypex',
                    type       => 'txt',
                    metadata   => { expected_sample => $sample, source_bam => $source_bam }
                );
                my $gtypex_path = $gtypex_file->path;
                my $cmd         = qq[$htscmd_exe gtcheck -g $genotypes_vcf $vcf_path > $gtypex_path];
                $self->dispatch([$cmd, $req, { output_files => [$gtypex_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            htscmd_gtcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'file of genotype concurrence scores calculated by htscmd gtcheck',
                metadata    => { expected_sample => 'name of expected sample', source_bam => 'input bam path' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs htscmd gtcheck to check sample identity in query vcf files against genotypes vcf files";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
