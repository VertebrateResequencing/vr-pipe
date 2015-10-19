
=head1 NAME

VRPipe::Steps::filter_unitig_pileup_calls - a step

=head1 DESCRIPTION

....

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::filter_unitig_pileup_calls with VRPipe::StepRole {
    method options_definition {
        return {
            k8_exe                  => VRPipe::StepOption->create(description => 'path to your k8 executable',                                                                                                             default_value => 'k8'),
            hapdip                  => VRPipe::StepOption->create(description => 'path to your hapdip.js javascript',                                                                                                      default_value => 'hapdip.js'),
            hapdip_deovlp_options   => VRPipe::StepOption->create(description => 'options to hapdip.js deovlp',                                                                                                            optional      => 1),
            hapdip_anno_options     => VRPipe::StepOption->create(description => 'options to hapdip.js anno',                                                                                                              optional      => 1),
            hapdip_filter_options   => VRPipe::StepOption->create(description => 'options to hapdip.js filter',                                                                                                            optional      => 1, default_value => '-q3'),
            post_hapdip_filter_cmds => VRPipe::StepOption->create(description => 'after running hapdip filter, option to pipe output vcf through another command (must contain either $output_vcf or $output_bcf string)', default_value => 'bgzip -c > $output_vcf'),
            bcftools_exe            => VRPipe::StepOption->create(description => 'path to bcftools executable (Will also be used to replace $bcftools placeholder in post_vep_options)',                                   default_value => 'bcftools'),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => '1 or more vcf or bcf files to filter'),
        };
    }
    
    method body_sub {
        return sub {
            my $self                  = shift;
            my $options               = $self->options;
            my $k8                    = $options->{k8_exe};
            my $hapdip                = $options->{hapdip};
            my $hapdip_deovlp_options = $options->{hapdip_deovlp_options};
            my $hapdip_anno_options   = $options->{hapdip_anno_options};
            my $hapdip_filter_options = $options->{hapdip_filter_options};
            my $post_filter_cmds      = $options->{post_hapdip_filter_cmds};
            my $bcftools              = $options->{bcftools_exe};
            
            if ($post_filter_cmds !~ /\$output_vcf/ && $post_filter_cmds !~ /\$output_bcf/) {
                $self->throw("post_hapdip_filter_cmds must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'hapdip.js',
                    version => VRPipe::StepCmdSummary->determine_version(qq[$k8 $hapdip], '^Version\s+(.+)$'),
                    summary => "k8 hapdip.js deovlp $hapdip_deovlp_options \$input_vcf | k8 hapdip.js anno $hapdip_anno_options | gzip -1 > tmp.vcf.gz && k8 hapdip.js filter $hapdip_filter_options tmp.vcf.gz | $post_filter_cmds"
                )
            );
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            
            my $suffix     = $post_filter_cmds =~ /\$output_bcf/ ? "bcf" : "vcf.gz";
            my $type       = $post_filter_cmds =~ /\$output_bcf/ ? "bcf" : "vcf";
            my $idx_suffix = $post_filter_cmds =~ /\$output_bcf/ ? "csi" : "tbi";
            foreach my $input_file (@{ $self->inputs->{vcf_files} }) {
                my $input_meta = $input_file->metadata;
                my $basename   = $input_file->basename;
                $basename =~ s/(\.vcf\.gz$|\.bcf$)//;
                my $temp_path    = $self->output_file(basename   => "$basename.tmp.vcf.gz",            type     => 'vcf',                           temporary => 1)->path;
                my $output_file  = $self->output_file(output_key => 'hapdip_filtered_vcf_files',       basename => "$basename.$suffix",             type      => $type, metadata => $input_meta);
                my $output_index = $self->output_file(output_key => 'hapdip_filtered_vcf_index_files', basename => "$basename.$suffix.$idx_suffix", type      => $idx_suffix, metadata => $input_meta);
                my $input_path  = $input_file->path;
                my $output_path = $output_file->path;
                my $post_filter = $post_filter_cmds;
                $post_filter =~ s/\$output_(b|v)cf/$output_path/g;
                $post_filter =~ s/\$bcftools/$bcftools/g;
                
                my $cmd = "$k8 $hapdip deovlp $hapdip_deovlp_options $input_path | $k8 $hapdip anno $hapdip_anno_options | gzip -1 > $temp_path && ";
                $cmd .= "$k8 $hapdip filter $hapdip_filter_options $temp_path | $post_filter && ";
                my $index = $type eq 'vcf' ? 'index -t' : 'index';
                $cmd .= "$bcftools $index $output_path";
                
                $self->dispatch([$cmd, $req, { output_files => [$output_file, $output_index] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            hapdip_filtered_vcf_files       => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => 'output .bcf or .vcf files'),
            hapdip_filtered_vcf_index_files => VRPipe::StepIODefinition->create(type => 'idx', max_files => -1, description => 'output CSI index for VCF/BCF files')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Filter fermikit unitig pileup calls using hapdip.js.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
