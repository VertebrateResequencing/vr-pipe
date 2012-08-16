
=head1 NAME

VRPipe::Steps::vcf_concat - a step

=head1 DESCRIPTION

Runs vcf-concat against an input set of VCFs, each containing a sequence number
as metadata for sorting purposes, generating one concatenated VCF

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

class VRPipe::Steps::vcf_concat with VRPipe::StepRole {
    method options_definition {
        return {
            vcf_concat_exe        => VRPipe::StepOption->create(description => 'path to vcf-concat executable',                                   optional => 1, default_value => 'vcf-concat'),
            vcf_concat_sites_only => VRPipe::StepOption->create(description => 'do not output genotype information to the concatenated vcf file', optional => 1, default_value => 0)
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'vcf files to concat',
                metadata    => { seq_no => 'a sequence number assigned by the split for reassembly in correct order' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self           = shift;
            my $options        = $self->options;
            my $vcf_concat_exe = $options->{vcf_concat_exe};
            my $sites_only     = $options->{vcf_concat_sites_only};
            
            # create temporary fofn of files to merge
            my $merge_list = $self->output_file(basename => "merge_list.txt", type => 'txt', temporary => 1);
            my @sorted_vcf_files = sort { $a->metadata->{seq_no} <=> $b->metadata->{seq_no} } @{ $self->inputs->{vcf_files} };
            $merge_list->create_fofn(\@sorted_vcf_files);
            
            # define output file
            my $concat_meta = $self->common_metadata($self->inputs->{vcf_files});
            my $concat_vcf = $self->output_file(output_key => 'concat_vcf', basename => "merged.vcf.gz", type => 'vcf', metadata => $concat_meta);
            
            # run command
            my $merge_list_path = $merge_list->path;
            my $concat_vcf_path = $concat_vcf->path;
            my $cut             = $sites_only ? ' | cut -f 1-8' : '';
            my $cmd             = qq[$vcf_concat_exe -f $merge_list_path$cut | bgzip -c > $concat_vcf_path];
            my $req             = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch([$cmd, $req, { output_files => [$concat_vcf, $merge_list] }]);
        };
    }
    
    method outputs_definition {
        return { concat_vcf => VRPipe::StepIODefinition->create(type => 'vcf', max_files => 1, description => 'a concatenated .vcf.gz file for each set of input vcfs') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run vcf-concat against input set of vcfs, each containing a sequence number as metadata, generating one concatenated vcf";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
