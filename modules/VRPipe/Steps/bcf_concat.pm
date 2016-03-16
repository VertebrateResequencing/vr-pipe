
=head1 NAME

VRPipe::Steps::bcf_concat - a step

=head1 DESCRIPTION

Runs vcf-concat against an input set of VCFs, each containing a sequence number
as metadata for sorting purposes, generating one concatenated VCF

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2016 Genome Research Limited.

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

class VRPipe::Steps::bcf_concat extends VRPipe::Steps::bcftools {
    around options_definition {
        return {
            %{ $self->$orig },
            bcftools_concat_options => VRPipe::StepOption->create(description => 'options to bcftools concat excluding -f; must include either $output_vcf or $output_bcf placeholder to define output file', default_value => '-Oz -o $output_vcf'),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'var',
                max_files   => -1,
                description => 'input VCF or BCF files to concat',
                metadata    => { seq_no => 'a sequence number assigned by the split for reassembly in correct order' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self           = shift;
            my $options        = $self->options;
            my $bcftools       = $options->{bcftools_exe};
            my $concat_options = $options->{bcftools_concat_options};
            
            if ($concat_options !~ /\$output_vcf/ && $concat_options !~ /\$output_bcf/) {
                $self->throw("concat_options must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            
            my $version = VRPipe::StepCmdSummary->determine_version($bcftools, '^Version: (.+)$');
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => $version,
                    summary => qq[bcftools concat -f \$concat.list $concat_options"],
                )
            );
            
            # create temporary fofn of files to merge
            my $concat_list = $self->output_file(basename => "concat_list.txt", type => 'txt', temporary => 1);
            my @sorted_vcf_files = sort { $a->metadata->{seq_no} <=> $b->metadata->{seq_no} } @{ $self->inputs->{vcf_files} };
            $concat_list->create_fofn(\@sorted_vcf_files);
            
            # define output file
            my $concat_meta = $self->common_metadata($self->inputs->{vcf_files});
            my $type        = $concat_options =~ /\$output_vcf/ ? 'vcf' : 'bcf';
            my $suffix      = $type eq 'vcf' ? 'vcf.gz' : 'bcf';
            my $idx         = $type eq 'vcf' ? 'tbi' : 'csi';
            my $basename    = "concat.$suffix";
            my $concat_vcf  = $self->output_file(output_key => 'concat_variant_file', basename => $basename, type => $type, metadata => $concat_meta);
            my $concat_idx  = $self->output_file(output_key => 'concat_variant_index_file', basename => "$basename.$idx", type => $idx, metadata => $concat_meta);
            
            # requirements
            my ($cpus) = $concat_options =~ m/--threads\s*(\d+)/;
            my $req = $self->new_requirements(memory => 1000, time => 1, $cpus ? (cpus => $cpus) : ());
            
            # run command
            my $concat_list_path = $concat_list->path;
            my $concat_vcf_path  = $concat_vcf->path;
            $concat_options =~ s/\$bcftools/$bcftools/g;
            $concat_options =~ s/\$output_[vb]cf/$concat_vcf_path/g;
            my $cmd = qq[$bcftools concat -f $concat_list_path $concat_options];
            my $idx_opt = $type eq 'vcf' ? '-ft' : '-f';
            $cmd = "($cmd) && $bcftools index $idx_opt $concat_vcf_path";
            $self->dispatch([$cmd, $req, { output_files => [$concat_vcf, $concat_idx] }]);
        };
    }
    
    method outputs_definition {
        return {
            concat_variant_file       => VRPipe::StepIODefinition->create(type => 'var', max_files => 1, description => 'a concatenated .vcf.gz file for each set of input vcfs'),
            concat_variant_index_file => VRPipe::StepIODefinition->create(type => 'idx', max_files => 1, description => 'a concatenated .vcf.gz file for each set of input vcfs'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run bcftools concat on an input set of VCFs or BCFs, generating one concatenated VCF or BCF file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
