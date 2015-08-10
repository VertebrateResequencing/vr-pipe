
=head1 NAME

VRPipe::Steps::vep_annotate_with_genome_chunking - a step

=head1 DESCRIPTION

Runs Ensembl's variant effect predictor to annotate input vcf files in chunks.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

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

class VRPipe::Steps::vep_annotate_with_genome_chunking extends VRPipe::Steps::vep_annotate with VRPipe::StepGenomeChunkingRole {
    method inputs_definition {
        return {
            vcf_file => VRPipe::StepIODefinition->create(type => 'var', max_files => 1, description => 'one vcf or bcf file to annotate'),
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $vep_exe           = $options->{vep_exe};
            my $vep_opts          = $options->{vep_options};
            my $bcftools_exe      = $options->{bcftools_exe};
            my $post_vep_cmds     = $options->{post_vep_options};
            my $ensembl_api_paths = $options->{ensembl_api_paths};
            if ($vep_opts =~ /-[i,o] /) {
                $self->throw("vep_options should not include the -i or -o option");
            }
            if (defined $post_vep_cmds && $post_vep_cmds !~ /\$output_vcf/ && $post_vep_cmds !~ /\$output_bcf/) {
                $self->throw("post_vep_options must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'variant_effect_predictor.pl',
                    version => VRPipe::StepCmdSummary->determine_version($vep_exe, '^version (.+)$'),
                    summary => "perl variant_effect_predictor.pl $vep_opts -i \$input -o \$output"
                )
            );
            my $suffix     = $post_vep_cmds =~ /\$output_bcf/ ? 'bcf' : 'vcf.gz';
            my $input_file = $self->inputs->{vcf_file}[0];
            my $input_path = $input_file->path;
            my $input_meta = $input_file->metadata;
            my $chunks     = $self->chunks();
            foreach my $chunk (@$chunks) {
                next if (exists $$input_meta{chrom} && $$input_meta{chrom} ne $chunk->{chrom});
                my $chrom          = $chunk->{chrom};
                my $from           = $chunk->{from};
                my $to             = $chunk->{to};
                my @chunk          = ($chrom, $from, $to);
                my $chunk_meta     = { %$input_meta, %$chunk };
                my $chunk_basename = "${chrom}_${from}-${to}.vep_annot.$suffix";
                #$self->output_file(basename => "${chunk_basename}_temp.vcf", type => 'vcf', temporary => 1);
                my $output_file  = $self->output_file(output_key => 'vep_annot_vcf_file',        basename => $chunk_basename,       type => 'var', metadata => $chunk_meta);
                my $output_index = $self->output_file(output_key => 'vep_annot_vcf_index_files', basename => "$chunk_basename.csi", type => 'idx', metadata => $chunk_meta);
                my $output_path  = $output_file->path;
                my $post_vep     = $post_vep_cmds;
                $post_vep =~ s/\$output_(b|v)cf/$output_path/;
                my $req = $self->new_requirements(memory => 5000, time => 1);
                my $this_cmd = "use VRPipe::Steps::vep_annotate; VRPipe::Steps::vep_annotate->vep_annotate_and_check(input => q[$input_path], chunk => [qw(@chunk)], output => q[$output_path], bcftools => q[$bcftools_exe], vep => q[$vep_exe], vep_opts => q[$vep_opts], api_paths => q[$ensembl_api_paths], post_vep => q[$post_vep]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$output_file, $output_index] });
            }
        };
    }
    
    method outputs_definition {
        return {
            vep_annot_vcf_file        => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => 'output .bcf or .vcf chunks'),
            vep_annot_vcf_index_files => VRPipe::StepIODefinition->create(type => 'idx', max_files => -1, description => 'output CSI index for VCF/BCF chunks')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run Ensembl VEP to annotate one or more vcf or bcf files. Split annotating into multiple jobs across the genome generating one vcf or bcf file per chunk.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
