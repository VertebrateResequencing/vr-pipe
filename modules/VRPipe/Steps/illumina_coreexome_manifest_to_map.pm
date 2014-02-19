
=head1 NAME

VRPipe::Steps::illumina_coreexome_manifest_to_map - a step

=head1 DESCRIPTION

This step processes an Illumina HumanCoreExome BeadChip manifest file to
generate a mapping file suitable for use with the genome_studio_fcr_to_vcf
step. (Both these steps use the same fcr-to-vcf exe found in our vr-codebase
git repository.)

This step requires you have the standard unix tool 'sort' installed, as well as
'bgzip' from the samtools repository.

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

class VRPipe::Steps::illumina_coreexome_manifest_to_map with VRPipe::StepRole {
    use DateTime;
    
    method options_definition {
        return {
            coreexome_manifest => VRPipe::StepOption->create(description => 'Illumina HumanCoreExome BeadChip manifest csv', optional => 0),
            fcr_to_vcf_exe     => VRPipe::StepOption->create(
                description   => 'path to the fcr-to-vcf executable',
                optional      => 1,
                default_value => 'fcr-to-vcf'
            ),
            bgzip_exe => VRPipe::StepOption->create(
                description   => 'path to the bgzip executable',
                optional      => 1,
                default_value => 'bgzip'
            ),
            sort_exe => VRPipe::StepOption->create(
                description   => 'path to the sort executable',
                optional      => 1,
                default_value => 'sort'
            )
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method outputs_definition {
        return {
            map_file => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'tabbed gz file'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $manifest_file = file($options->{coreexome_manifest});
            $self->throw("coreexome_manifest file $manifest_file does not exist") unless -s $manifest_file;
            
            my $fcr_to_vcf = $options->{fcr_to_vcf_exe};
            my $sort       = $options->{sort_exe};
            my $bgzip      = $options->{bgzip_exe};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'fcr-to-vcf',
                    version => 0,
                    summary => "cat $manifest_file | $fcr_to_vcf -c | $sort -t, -k10,10d -k11,11n | $fcr_to_vcf -M | $bgzip -c > \$output_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $ofile = $self->output_file(output_key => 'map_file', output_dir => $manifest_file->dir, basename => $manifest_file->basename . '.map.tab.gz', type => 'bin', metadata => { manifest_file => $manifest_file })->path;
            $self->dispatch(["cat $manifest_file | $fcr_to_vcf -c | $sort -t, -k10,10d -k11,11n | $fcr_to_vcf -M | $bgzip -c > $ofile", $req, { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method description {
        return "Convert an Illumina HumanCoreExome BeadChip manifest file to a mapping suitable for helping convert fcr files to vcf.";
    }
}

1;
