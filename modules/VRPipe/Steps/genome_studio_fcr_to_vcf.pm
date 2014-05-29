
=head1 NAME

VRPipe::Steps::genome_studio_fcr_to_vcf - a step

=head1 DESCRIPTION

This step converts the data in a single-sample fcr file (eg. produced by the
split_genome_studio_genotype_files step) into an accurate VCF with correct SNP
positions and strand orientation.

It uses the fcr-to-vcf exe from the vr-codebase git repository.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2014 Genome Research Limited.

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

class VRPipe::Steps::genome_studio_fcr_to_vcf with VRPipe::StepRole {
    use DateTime;
    
    method options_definition {
        return {
            vcf_sample_from_metadata => VRPipe::StepOption->create(
                description   => 'if the id in the second column of the fcr matches metadata with key x, but you want the sample id in the VCF to have metadata from key y, provide x:y; separate multiple y keys with + symbols - values will be joined with underscores',
                optional      => 1,
                default_value => 'infinium_sample:public_name+sample'
            ),
            fcr_to_vcf_exe => VRPipe::StepOption->create(
                description   => 'path to the fcr-to-vcf executable',
                optional      => 1,
                default_value => 'fcr-to-vcf'
            ),
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to the bcftools executable (used by fcr-to-vcf)',
                optional      => 1,
                default_value => 'bcftools'
            )
        };
    }
    
    method inputs_definition {
        return {
            map_file => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'mapping file produced by the illumina_coreexome_manifest_to_map step'
            ),
            fcr_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'single-sample fcr file with sample metadata',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line' }
            ),
        };
    }
    
    method outputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'VCF file',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line' }
            ),
            tbi_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'tbi file',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self       = shift;
            my $options    = $self->options;
            my $fcr_to_vcf = $options->{fcr_to_vcf_exe};
            my $bcftools   = $options->{bcftools_exe};
            my $sfm        = $options->{vcf_sample_from_metadata};
            my ($src_key, @dst_keys);
            my $summary_s = '';
            if ($sfm) {
                ($src_key, my $dst_key) = split(':', $sfm);
                @dst_keys = split(/\+/, $dst_key);
                $summary_s = " -s $sfm";
            }
            else {
                $summary_s = " -m";
            }
            
            my ($map_file) = @{ $self->inputs->{map_file} };
            $map_file = $map_file->path;
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'fcr-to-vcf',
                    version => 0,
                    summary => "cat \$fcr_file | $fcr_to_vcf -b $bcftools -a \$map_file$summary_s -o \$outdir"
                )
            );
            
            foreach my $fcr_file (@{ $self->inputs->{fcr_files} }) {
                my $meta = $fcr_file->metadata;
                
                my $outdir = $meta->{sample};
                my $basename;
                my $s = '';
                if ($src_key && @dst_keys && defined $meta->{$src_key} && defined $meta->{ $dst_keys[0] }) {
                    my $dst_vals = join('_', map { $meta->{$_} || 'undef' } @dst_keys);
                    $s        = " -s $meta->{$src_key}:$dst_vals";
                    $outdir   = $dst_vals;
                    $basename = file($outdir, "$outdir.vcf.gz")->stringify;
                }
                else {
                    $s        = ' -m';
                    $basename = "$outdir.vcf.gz";
                }
                
                #$self->output_file(basename => $basename, type => 'vcf', temporary => 1);
                
                my $vcf_file_path = $self->output_file(output_key => 'vcf_files', basename => $basename, type => 'vcf', metadata => $meta)->path;
                $self->output_file(output_key => 'tbi_files', basename => $basename . '.tbi', type => 'tbi');
                my $fcr_file_path = $fcr_file->path;
                
                $self->dispatch(["cat $fcr_file_path | $fcr_to_vcf -b $bcftools -a $map_file$s -o $outdir", $req]);
            }
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method description {
        return "Convert single-sample genotype data text files (FCR files) to sorted compressed VCF suitable for calling with.";
    }
}

1;
