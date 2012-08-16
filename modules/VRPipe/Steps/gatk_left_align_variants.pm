
=head1 NAME

VRPipe::Steps::gatk_left_align_variants - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

#Example generic command for UnifiedGenotyper GATK v1.3
#java -Xmx2g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T LeftAlignVariants \
#   --variant input.vcf \
#   -o output.vcf

class VRPipe::Steps::gatk_left_align_variants extends VRPipe::Steps::gatk {
    around options_definition {
        return {
            %{ $self->$orig },
            left_align_variants_options => VRPipe::StepOption->create(description => 'Any addional general GATK options to pass the LeftAlignVariants', optional => 1),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'input vcf files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $reference_fasta   = $options->{reference_fasta};
            my $leftalign_options = $options->{left_align_variants_options};
            
            my $req = $self->new_requirements(memory => 1200, time => 1);
            my $jvm_args = $self->jvm_args($req->memory);
            
            my $idx = 0;
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $basename = $vcf->basename;
                $basename =~ s/vcf(\.gz)?/aln.$idx.vcf.gz/;
                my $vcf_out     = $self->output_file(output_key => 'left_aligned_vcf_files', basename => $basename, type => 'vcf');
                my $input_path  = $vcf->path;
                my $output_path = $vcf_out->path;
                my $cmd         = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T LeftAlignVariants -R $reference_fasta --variant $input_path -o $output_path $leftalign_options];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_left_align_variants', 'left_align_variants', [$cmd, $req, { output_files => [$vcf_out] }]);
                ++$idx;
            }
        };
    }
    
    method outputs_definition {
        return { left_aligned_vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'output vcf files') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs GATK LeftAlignVariants to left-align indels in VCF files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method left_align_variants (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /-variant (\S+) \-o (\S+)/;
        
        my $input_file = VRPipe::File->get(path => $input_path);
        my $input_recs = $input_file->num_records;
        $input_file->disconnect;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_recs = $output_file->num_records;
        
        unless ($output_recs == $input_recs) {
            $output_file->unlink;
            $self->throw("Output VCF has different number of data lines from input (input $input_recs, output $output_recs)");
        }
        else {
            return 1;
        }
    }
}

1;
