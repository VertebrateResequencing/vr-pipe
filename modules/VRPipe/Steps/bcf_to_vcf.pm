=head1 NAME

VRPipe::Steps::bcf_to_vcf - a step

=head1 DESCRIPTION

*** more documentation to come

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

class VRPipe::Steps::bcf_to_vcf with VRPipe::StepRole {
    method options_definition {
        return { bcftools_exe => VRPipe::StepOption->create(description => 'path to bcftools executable', optional => 1, default_value => 'bcftools'),
                 bcftools_view_options => VRPipe::StepOption->create(description => 'bcftools view options', optional => 1, default_value => '-p 0.99 -vcgN') };
    }
    method inputs_definition {
        return { bcf_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => '1 or more bcf files to convert to compressed vcf'),
                 samples_files => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0, max_files => -1, description => 'Optional samples file for restricting samples to call on and/or defining ploidy for each sample'),
                 # sites_file => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0, max_files => 1, description => 'Optional sites file for calling only at the given sites'),
             };
    }

    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $bcftools = $options->{bcftools_exe};
            my $view_opts = $options->{bcftools_view_options};
            
            my %samples;
            foreach my $sample_file (@{$self->inputs->{samples_files}}) {
                my $bcf = $sample_file->metadata->{source_bcf} || next;
                $samples{$bcf} = $sample_file->path;
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bcf (@{$self->inputs->{bcf_files}}) {
                my $bcf_meta = $bcf->metadata;
                my $bcf_path = $bcf->path;
                my $basename = $bcf->basename;
                $basename =~ s/bcf$/vcf.gz/;
                my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => $basename, type => 'vcf', metadata => $bcf_meta);
                my $vcf_path = $vcf_file->path;
                
                my $sample_opts = '';
                if (exists $samples{$bcf_path}) {
                    $sample_opts = " -s $samples{$bcf_path}";
                }
                my $cmd = qq[$bcftools view $view_opts$sample_opts $bcf_path | bgzip -c > $vcf_path];
                $self->dispatch([$cmd, $req, {output_files => [$vcf_file]}]); 
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'bcftools', 
                                   version => VRPipe::StepCmdSummary->determine_version($bcftools, '^Version: (.+)$'), 
                                   summary => "bcftools view $view_opts \$bcf_file | bgzip -c > \$vcf_file"));
        };
    }
    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'a .vcf.gz file for each input bcf file') };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Run bcftools view option to generate one compressed vcf file per input bcf";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
