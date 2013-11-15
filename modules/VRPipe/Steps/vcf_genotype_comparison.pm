
=head1 NAME

VRPipe::Steps::vcf_genotype_comparison - a step

=head1 DESCRIPTION

Compare the genotypes of the samples in a VCF using bcftools gtcheck. Adds
'genotype_maximum_deviation' metadata to the VCF file indicating if all its
samples were similar to each other.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::vcf_genotype_comparison with VRPipe::StepRole {
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to the bcftools executable',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_gtcheck_options => VRPipe::StepOption->create(
                description => 'options to bcftools gtcheck (-g option is not valid for this step)',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'vcf files; genotype comparison will be done on each independently'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self         = shift;
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            my $gtcheck_opts = $options->{bcftools_gtcheck_options};
            if ($gtcheck_opts =~ /gtcheck| -g| --genotypes/) {
                $self->throw("bcftools_gtcheck_options should not include the gtcheck subcommand or the -g option");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => VRPipe::StepCmdSummary->determine_version($bcftools_exe, '^Version: (.+)$'),
                    summary => "bcftools gtcheck $gtcheck_opts \$vcf_file > \$gtcheck_file.gtypex"
                )
            );
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $vcf_path = $vcf->path;
                my $meta     = $vcf->metadata;
                
                my $gtypex_file = $self->output_file(
                    output_key => 'bcftools_gtcheck_files',
                    basename   => $vcf->basename . '.gtypex',
                    type       => 'txt',
                    metadata   => $vcf->metadata
                );
                my $gtypex_path = $gtypex_file->path;
                my $cmd         = qq[$bcftools_exe gtcheck $gtcheck_opts $vcf_path > $gtypex_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_genotype_comparison', 'compare_genotypes', [$cmd, $req, { output_files => [$gtypex_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bcftools_gtcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'file of genotype concurrence scores calculated by bcftools gtcheck',
                metadata    => { genotype_maximum_deviation => "maximum deviation in genotype and the sample causing that deviation" }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Compare the genotypes of the samples in a VCF to each other to confirm they come from the same individual using bcftools gtcheck";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method compare_genotypes (ClassName|Object $self: Str $cmd_line) {
        my ($vcf_path, $output_path) = $cmd_line =~ /(\S+) > (\S+)$/;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $vcf_file    = VRPipe::File->get(path => $vcf_path);
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        
        my $fh = $output_file->openr;
        while (<$fh>) {
            next unless /^MD\s+(\S+)\s+(\S+)/;
            my $md     = $1;
            my $sample = $2;
            
            foreach my $file ($vcf_file, $output_file) {
                $file->add_metadata({ genotype_maximum_deviation => "$md:$sample" });
            }
            
            last;
        }
        $output_file->close;
    }
}

1;
