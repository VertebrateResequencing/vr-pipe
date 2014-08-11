
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
            
            if ($bcfopts =~ /-O\s*[zv]/) {
                # we're hard-coded for bcf output, so we really don't support
                # these vcf options
                $self->throw("-O z or v are not supported - this step can only output bcfs");
            }
            
            my $merged_basename = 'merged.bcf';
            my $merged_meta     = $self->combined_metadata($self->inputs->{vcf_files});
            $self->_merge($bcftools_exe, $bcfopts, $self->inputs->{vcf_files}, $merged_basename, 'merged_bcf', 'bcf', $merged_meta, 'bcf_index');
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
