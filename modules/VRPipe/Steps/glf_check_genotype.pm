
=head1 NAME

VRPipe::Steps::glf_check_genotype - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Steps::glf_check_genotype with VRPipe::StepRole {
    method options_definition {
        return {
            hapmap2bin_sample_genotypes_file => VRPipe::StepOption->create(description => 'absolute path to binary file of sample genotypes, produced by hapmap2bin'),
            expected_sample_metadata_key     => VRPipe::StepOption->create(
                description    => 'The key of the metadata on the bam/bcf file that holds the expected sample name to compare against in the hapmap2bin_sample_genotypes_file file',
                optional       => 1,
                default_value  => 'individual',
                allowed_values => ['individual', 'sample']
            ),
            glf_exe => VRPipe::StepOption->create(
                description   => 'path to the glf executable',
                optional      => 1,
                default_value => 'glf'
            )
        };
    }
    
    method inputs_definition {
        return {
            bcf_files => VRPipe::StepIODefinition->create(
                type        => 'bcf',
                max_files   => -1,
                description => 'bcf files for genotyping',
                metadata    => {
                    source_bam => 'input bam path',
                    individual => 'expected individual name',
                    sample     => 'expected sample name'
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self           = shift;
            my $options        = $self->options;
            my $gtype_snps_bin = Path::Class::File->new($options->{hapmap2bin_sample_genotypes_file});
            $self->throw("hapmap2bin_sample_genotypes_file must be an absolute path") unless $gtype_snps_bin->is_absolute;
            my $glf_exe      = $options->{glf_exe};
            my $expected_key = $options->{expected_sample_metadata_key};
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $bcf (@{ $self->inputs->{bcf_files} }) {
                my $bcf_path    = $bcf->path;
                my $meta        = $bcf->metadata;
                my $sample      = $meta->{$expected_key};
                my $source_bam  = $meta->{source_bam};
                my $gtypex_file = $self->output_file(
                    output_key => 'gtypex_files_with_metadata',
                    basename   => $bcf->basename . '.gtypex',
                    type       => 'txt',
                    metadata   => { expected_sample => $sample, source_bam => $source_bam }
                );
                my $gtypex_path = $gtypex_file->path;
                my $cmd         = qq[$glf_exe checkGenotype -s - $gtype_snps_bin $bcf_path > $gtypex_path];
                $self->dispatch([$cmd, $req, { output_files => [$gtypex_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            gtypex_files_with_metadata => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'file of likelihood scores calculated by glf',
                metadata    => { expected_sample => 'name of expected sample', source_bam => 'input bam path' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Produces gtypex files of genotype likelihood scores for sample genotypes which are calculated by running glf checkGenotype against a snp binary file (of all samples) on bcf files";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
