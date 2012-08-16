
=head1 NAME

VRPipe::Steps::gatk_variant_annotator - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

# java -Xmx2g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T VariantAnnotator \
#   -I input.bam \
#   -o output.vcf \
#   -A DepthOfCoverage \
#   --variant input.vcf \
#   -L input.vcf \
#   --dbsnp dbsnp.vcf

class VRPipe::Steps::gatk_variant_annotator extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{ $self->$orig }, variant_annotator_options => VRPipe::StepOption->create(description => 'Options for GATK VariantAnnotator.'), };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'one or more vcf files to be annotated'),
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', min_files => 0,  max_files   => -1, description => 'one or more bam files')
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $reference_fasta = $options->{reference_fasta};
            my $annotation_opts = $options->{variant_annotator_options};
            if ($annotation_opts =~ /VariantAnnotator|-variant/) {
                $self->throw("variant_annotator_options should not include the VariantAnnotator task command or the -variant option");
            }
            
            my $req = $self->new_requirements(memory => 1200, time => 1);
            my $jvm_args = $self->jvm_args($req->memory);
            
            my $bams_list;
            if ($self->inputs->{bam_files}) {
                $bams_list = $self->output_file(basename => "bams.list", type => 'txt', temporary => 1);
                $bams_list->create_fofn($self->inputs->{bam_files});
            }
            
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $inputs = '--variant ' . $vcf->path;
                if ($bams_list) {
                    $inputs .= ' -I ' . $bams_list->path;
                }
                my $basename = $vcf->basename;
                $basename =~ s/vcf(.gz)?$/anno.vcf.gz/;
                
                my $anno_vcf_file = $self->output_file(output_key => 'annotated_vcf_files', basename => $basename, type => 'vcf', metadata => $vcf->metadata);
                my $anno_vcf_path = $anno_vcf_file->path;
                
                my $cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T VariantAnnotator -R $reference_fasta $inputs -o $anno_vcf_path $annotation_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_variant_annotator', 'annotate_variants', [$cmd, $req, { output_files => [$anno_vcf_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { annotated_vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'An annotated VCF files') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run GATK VariantAnnotator to annotate multiple VCF files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method annotate_variants (ClassName|Object $self: Str $cmd_line) {
        my ($input_path)  = $cmd_line =~ /--variant (\S+)/;
        my ($output_path) = $cmd_line =~ /-o (\S+)/;
        my $input_file  = VRPipe::File->get(path => $input_path);
        my $output_file = VRPipe::File->get(path => $output_path);
        
        $output_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        $output_file->update_stats_from_disc(retries => 3);
        
        my $input_records  = $input_file->num_records;
        my $output_records = $output_file->num_records;
        
        unless ($input_records == $output_records) {
            $output_file->unlink;
            $self->throw("Number of records in output VCF, $output_records, not equal to number of records in input VCF, $input_records");
        }
        else {
            return 1;
        }
    }
}

1;
