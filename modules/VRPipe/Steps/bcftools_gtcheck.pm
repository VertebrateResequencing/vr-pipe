
=head1 NAME

VRPipe::Steps::bcftools_gtcheck - a step

=head1 DESCRIPTION

Runs bcftools gtcheck to check sample identities in query vcf files against
genotypes in vcf file.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bcftools_gtcheck extends VRPipe::Steps::bcftools {
    around options_definition {
        return {
            %{ $self->$orig },
            bcftools_gtcheck_options => VRPipe::StepOption->create(
                description => 'options to bcftools gtcheck (excluding -g)',
                optional    => 1
            ),
            expected_sample_from_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key that expected_sample should be taken from',
                optional      => 1,
                default_value => 'sample'
            ),
            genotypes_file => VRPipe::StepOption->create(
                description => 'absolute path to genotypes vcf or bcf file containing genotyping information on all possible samples; optional if you supply this via the datasource input',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'vcf files; genotype check will be done on each independently'
            ),
            genotypes_bcf => VRPipe::StepIODefinition->create(
                type        => 'bcf',
                min_files   => 0,
                max_files   => 1,
                description => 'bcf file containing genotyping information on all possible samples; optional if genotypes_file option is supplied'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self         = shift;
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            my $gtcheck_opts = $options->{bcftools_gtcheck_options};
            my $expected_key = $options->{expected_sample_from_metadata_key};
            $gtcheck_opts ||= '';
            if ($gtcheck_opts =~ /gtcheck| -g| --genotypes/) {
                $self->throw("bcftools_gtcheck_options should not include the gtcheck subcommand or the -g option");
            }
            
            my $genotypes_file_path;
            if ($options->{genotypes_file}) {
                $genotypes_file_path = file($options->{genotypes_file});
                $self->throw("genotypes_file must be an absolute path") unless $genotypes_file_path->is_absolute;
            }
            else {
                my ($genotypes_bcf) = @{ $self->inputs->{genotypes_bcf} || [] };
                $genotypes_file_path = $genotypes_bcf->path if $genotypes_bcf;
            }
            $genotypes_file_path || $self->throw("A genotypes file is required, either via the genotypes_file option or as an input from the datasource");
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => $self->bcftools_version_string,
                    summary => "bcftools gtcheck -g $genotypes_file_path $gtcheck_opts \$vcf_file > \$gtcheck_file.gtypex"
                )
            );
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $vcf_path   = $vcf->path;
                my $meta       = $vcf->metadata;
                my $sample     = $meta->{$expected_key};
                my $source_bam = $meta->{source_bam};
                
                my $gtypex_file = $self->output_file(
                    output_key => 'bcftools_gtcheck_files',
                    basename   => $vcf->basename . '.gtypex',
                    type       => 'txt',
                    metadata   => { expected_sample => $sample, source_bam => $source_bam }
                );
                my $gtypex_path = $gtypex_file->path;
                my $cmd         = qq[$bcftools_exe gtcheck -g $genotypes_file_path $gtcheck_opts $vcf_path > $gtypex_path];
                $self->dispatch([$cmd, $req, { output_files => [$gtypex_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bcftools_gtcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'file of genotype concurrence scores calculated by bcftools gtcheck',
                metadata    => { expected_sample => 'name of expected sample', source_bam => 'input bam path' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Compare the genotype of a sample in a VCF file to those in a genotypes VCF file using bcftools gtcheck";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
