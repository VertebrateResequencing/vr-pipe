
=head1 NAME

VRPipe::Steps::vcf_merge_different_samples_control_aware - a step

=head1 DESCRIPTION

Like vcf_merge_different_samples, but aware of which input VCF is the
'control'. Useful for any situation where you have multiple samples from the
same individual and want to ensure that one of the samples (eg. the control
sample) comes first in the merged file.

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

class VRPipe::Steps::vcf_merge_different_samples_control_aware extends VRPipe::Steps::vcf_merge_different_samples {
    around options_definition {
        return {
            %{ $self->$orig },
            control_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key to check on the input files to see which one is the control, and the key used on the output file to identify the control',
                default_value => 'sample_control'
            ),
            control_metadata_regex => VRPipe::StepOption->create(
                description   => 'a (case-insensitive) regular expression that matches the value of the control_metadata_key when the input VCF is for the control sample',
                default_value => '1'
            ),
            control_metadata_sample_key => VRPipe::StepOption->create(
                description   => 'the metadata key to get the sample from the input control file, to apply as metadata on the output file; use + symbols to comprise it of multiple values concatenated with underscores',
                default_value => 'public_name+sample'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options       = $self->options;
            my $bcftools_exe  = $options->{bcftools_exe};
            my $bcfopts       = $options->{bcftools_options};
            my $control_key   = $options->{control_metadata_key};
            my $control_regex = $options->{control_metadata_regex};
            my $sample_key    = $options->{control_metadata_sample_key};
            my @sample_keys   = split(/\+/, $sample_key);
            
            my @input_set;
            my $control_sample;
            foreach my $vcf_file (sort { $a->path cmp $b->path } @{ $self->inputs->{vcf_files} }) {
                # there's only supposed to be one control sample, but in case
                # there's more than one we just pick the first alphabetically
                unless ($control_sample) {
                    my $control_value = $vcf_file->meta_value($control_key);
                    if ($control_value && $control_value =~ /$control_regex/i) {
                        unshift(@input_set, $vcf_file->path);
                        $control_sample = join('_', map { $vcf_file->meta_value($_) } @sample_keys);
                        next;
                    }
                }
                
                push @input_set, $vcf_file->path;
            }
            $control_sample || $self->throw("There was no control sample identified amongst the input vcf files (@input_set).");
            
            my $merged_meta = $self->common_metadata($self->inputs->{vcf_files});
            $merged_meta->{$control_key} = $control_sample;
            
            # the rest of this body_sub is identical to parent *** can we avoid
            # this code duplication?...
            
            my $merged_basename = 'merged.vcf.gz';
            my $merged_vcf      = $self->output_file(output_key => 'merged_vcf', basename => $merged_basename, type => 'vcf', metadata => $merged_meta);
            my $output_path     = $merged_vcf->path;
            
            if (@input_set == 1) {
                # merge doesn't work on 1 input file; just symlink the input to
                # output
                my $source = $self->inputs->{vcf_files}->[0];
                $source->symlink($merged_vcf);
                $merged_vcf->metadata($merged_meta);
                $merged_vcf->update;
                return 1;
            }
            
            my $this_cmd = qq[$bcftools_exe merge $bcfopts @input_set > $output_path];
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => VRPipe::StepCmdSummary->determine_version($bcftools_exe, '^Version: (.+)$'),
                    summary => "bcftools merge $bcfopts \@input_vcfs > \$output_path"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_merge_different_samples_control_aware', 'merge_vcf', [$this_cmd, $req, { output_files => [$merged_vcf] }]);
        };
    }
    
    method description {
        return "Merges compressed VCFs using bcftools merge which contain different samples to produce a single VCF containing all the input samples, with the first sample being the one from the VCF tagged with the right control metadata, and identifying that sample as a control in metadata on the output file";
    }
}

1;
