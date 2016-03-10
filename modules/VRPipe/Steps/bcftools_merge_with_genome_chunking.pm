
=head1 NAME

VRPipe::Steps::bcftools_merge_with_genome_chunking - a step

=head1 DESCRIPTION

...

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

class VRPipe::Steps::bcftools_merge_with_genome_chunking extends VRPipe::Steps::bcftools with VRPipe::StepGenomeChunkingRole {
    around options_definition {
        return {
            %{ $self->$orig },
            bcftools_merge_options => VRPipe::StepOption->create(
                description   => 'bcftools merge options exluding -l/--file-list; must include $output_vcf or $output_bcf to indicate if command should expect VCF or BCF output. e.g. "-Ou --gvcf | $bcftools call -mv -Ob $output_bcf"',
                optional      => 1,
                default_value => '-Oz $output_vcf'
            ),
            reference_fasta => VRPipe::StepOption->create(
                description => 'absolute path to genome reference fasta file to be optionally used to replace $reference_fasta in the bcftools_merge_options command line',
                optional    => 1,
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files       => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => '1 or more VCF or BCF files to merge'),
            vcf_index_files => VRPipe::StepIODefinition->create(type => 'idx', max_files => -1, description => '1 or more TBI or CSI index files for the input VCF/BCF files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self            = shift;
            my $inputs          = $self->inputs->{vcf_files};
            my $merged_metadata = $self->common_metadata($inputs);
            $merged_metadata = { %$merged_metadata, $self->element_meta };
            my $options = $self->handle_override_options($merged_metadata);
            
            my $bcftools        = $options->{bcftools_exe};
            my $merge_opts      = $options->{bcftools_merge_options};
            my $reference_fasta = $options->{reference_fasta};
            
            if ($merge_opts !~ /\$output_vcf/ && $merge_opts !~ /\$output_bcf/) {
                $self->throw("merge_opts must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            
            my $cmd_line = qq[\$bcftools merge -r \${CHROM}:\${FROM}-\${TO} -l \$file_list $merge_opts];
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => $self->bcftools_version_string,
                    summary => $cmd_line,
                )
            );
            $cmd_line =~ s/\$bcftools/$bcftools/g;
            $cmd_line =~ s/\$reference_fasta/$reference_fasta/g if ($reference_fasta);
            
            my $file_list = $self->output_file(basename => "merge.list", type => 'txt', temporary => 1);
            my $file_path = $file_list->path;
            $cmd_line =~ s/\$file_list/$file_path/g;
            $file_list->create_fofn($inputs);
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            
            my $chunks   = $self->chunks();
            my $type     = $cmd_line =~ /\$output_vcf/ ? 'vcf' : 'bcf';
            my $suffix   = $type eq 'vcf' ? 'vcf.gz' : 'bcf';
            my $idx      = $type eq 'vcf' ? 'tbi' : 'csi';
            my $basename = "merge.$suffix";
            foreach my $chunk (@$chunks) {
                next if (exists $$merged_metadata{chrom} && $$merged_metadata{chrom} ne $chunk->{chrom});
                my $chrom = $chunk->{chrom};
                my $from  = $chunk->{from};
                my $to    = $chunk->{to};
                
                my $chunk_basename = "${chrom}_${from}-${to}.$basename";
                
                my $this_cmd_line = $cmd_line;
                $this_cmd_line =~ s/\${CHROM}/$chrom/g;
                $this_cmd_line =~ s/\${FROM}/$from/g;
                $this_cmd_line =~ s/\${TO}/$to/g;
                
                my $chunk_meta = { %$merged_metadata, %$chunk };
                my $vcf_file  = $self->output_file(output_key => 'merged_variant_files',       basename => $chunk_basename,        type => $type, metadata => $chunk_meta);
                my $vcf_index = $self->output_file(output_key => 'merged_variant_index_files', basename => "$chunk_basename.$idx", type => $idx,  metadata => $chunk_meta);
                
                my $output_path = $vcf_file->path;
                $this_cmd_line =~ s/\$output_[bv]cf/$output_path/g;
                my $idx_opt = $type eq 'vcf' ? '-ft' : '-f';
                $this_cmd_line = "($this_cmd_line) && $bcftools index $idx_opt $output_path";
                $self->dispatch([$this_cmd_line, $req, { output_files => [$vcf_file, $vcf_index] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            merged_variant_files       => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => 'a .vcf.gz or bcf file for each input bcf file'),
            merged_variant_index_files => VRPipe::StepIODefinition->create(type => 'idx', max_files => -1, description => 'a .vcf.gz file for each input bcf file'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run bcftools merge.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
