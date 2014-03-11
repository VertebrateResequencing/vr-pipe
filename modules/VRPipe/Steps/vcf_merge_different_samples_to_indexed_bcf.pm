
=head1 NAME

VRPipe::Steps::vcf_merge_different_samples_to_indexed_bcf - a step

=head1 DESCRIPTION

Run bcftools merge to merge a set of VCFs each with a different sample,
creating a single output BCF that is indexed, with all the input samples.

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

class VRPipe::Steps::vcf_merge_different_samples_to_indexed_bcf extends VRPipe::Steps::vcf_merge_different_samples {
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools executable',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_options => VRPipe::StepOption->create(
                description   => 'Options for the bcftools merge command. Does not support the --print-header option; -O z or v are not supported.',
                optional      => 1,
                default_value => '-O b'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            my $bcfopts      = $options->{bcftools_options};
            
            my @input_set;
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                push @input_set, $vcf_file->path;
            }
            my $merged_basename = 'merged.bcf';
            
            my $merged_meta      = $self->common_metadata($self->inputs->{vcf_files});
            my $merged_bcf       = $self->output_file(output_key => 'merged_bcf', basename => $merged_basename, type => 'bcf', metadata => $merged_meta);
            my $output_path      = $merged_bcf->path;
            my $merged_bcf_index = $self->output_file(output_key => 'bcf_index', basename => $merged_basename . '.csi', type => 'bin', metadata => $merged_meta);
            
            my $this_cmd;
            if (@input_set == 1) {
                # merge doesn't work on 1 input file; just convert the vcf to
                # a bcf
                $this_cmd = qq[$bcftools_exe view -Sb @input_set > $output_path && $bcftools_exe index $output_path];
                
                $self->set_cmd_summary(
                    VRPipe::StepCmdSummary->create(
                        exe     => 'bcftools',
                        version => VRPipe::StepCmdSummary->determine_version($bcftools_exe, '^Version: (.+)$'),
                        summary => "bcftools view -Sb \@input_vcfs > \$output_path && bcftools index \$output_path"
                    )
                );
            }
            else {
                $this_cmd = qq[$bcftools_exe merge $bcfopts @input_set > $output_path && $bcftools_exe index $output_path];
                
                $self->set_cmd_summary(
                    VRPipe::StepCmdSummary->create(
                        exe     => 'bcftools',
                        version => VRPipe::StepCmdSummary->determine_version($bcftools_exe, '^Version: (.+)$'),
                        summary => "bcftools merge $bcfopts \@input_vcfs > \$output_path && bcftools index \$output_path"
                    )
                );
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_merge_different_samples_to_indexed_bcf', 'merge_vcf', [$this_cmd, $req, { output_files => [$merged_bcf, $merged_bcf_index] }]);
        };
    }
    
    method outputs_definition {
        return {
            merged_bcf => VRPipe::StepIODefinition->create(
                type        => 'bcf',
                description => 'a merged bcf file',
                max_files   => 1
            ),
            bcf_index => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'index of the merged bcf file',
                max_files   => 1
            )
        };
    }
    
    method description {
        return "Merges compressed VCFs using bcftools merge which contain different samples to produce a single indexed BCF containing data for all the input samples";
    }
}

1;
