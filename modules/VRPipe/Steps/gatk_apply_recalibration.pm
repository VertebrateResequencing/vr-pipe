
=head1 NAME

VRPipe::Steps::gatk_apply_recalibration - a step

=head1 DESCRIPTION

*** more documentation to come

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

#Example ApplyRecalibration command - GATK v1.3
#java -Xmx3g -jar GenomeAnalysisTK.jar \
#   -T ApplyRecalibration \
#   -R reference/human_g1k_v37.fasta \
#   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
#   --ts_filter_level 99.0 \
#   -tranchesFile path/to/output.tranches \
#   -recalFile path/to/output.recal \
#   -o path/to/output.recalibrated.filtered.vcf

class VRPipe::Steps::gatk_apply_recalibration extends VRPipe::Steps::gatk {
    around options_definition {
        return {
            %{ $self->$orig },
            apply_recalibration_options => VRPipe::StepOption->create(description => 'options for GATK ApplyRecalibration, excluding reference genome, mode, input and output files'),
            apply_recalibration_mode    => VRPipe::StepOption->create(description => 'mode', optional => 1, default_value => 'BOTH'),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files          => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'one or more tabixed vcf files processed by VariantRecalibrator'),
            recalibration_file => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1,  description => 'a recalibration table file in CSV format for each input vcf'),
            tranches_file      => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1,  description => 'a tranches file generated by VariantRecalibrator for each vcf'),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref              = $options->{reference_fasta};
            my $apply_recal_opts = $options->{apply_recalibration_options};
            if ($apply_recal_opts =~ /$ref|ApplyRecalibration|-mode|-input|-output/) {
                $self->throw("apply_recalibration_options should not include the reference, mode, input or output options or ApplyRecalibration task command");
            }
            my $mode = $options->{apply_recalibration_mode};
            
            my $req = $self->new_requirements(memory => 1200, time => 1);
            
            my $recal_file    = $self->inputs->{recalibration_file}[0];
            my $tranches_file = $self->inputs->{tranches_file}[0];
            
            my $idx = 0;
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $vcf_path = $vcf->path;
                my $vcf_meta = $vcf->metadata || {};
                my $basename = $vcf->basename;
                $basename =~ s/\.vcf(.gz)?$/.recal_$idx.vcf.gz/;
                
                my $recal_file_path    = $recal_file->path;
                my $tranches_file_path = $tranches_file->path;
                
                # If VCFs are split by chromosome, save time by only walking along the chromosome
                # Add interval set rule, so that if apply recalibration options contains -L option too (eg exome targets file)
                # then will work on the intersection
                my $this_apply_recal_opts = $apply_recal_opts;
                if (defined $vcf_meta->{chrom}) {
                    $this_apply_recal_opts .= " -L $$vcf_meta{chrom} --interval_set_rule INTERSECTION";
                    $basename = "$$vcf_meta{chrom}.recal_$idx.vcf.gz";
                }
                
                my $vcf_recal_file = $self->output_file(output_key => 'recalibrated_vcfs', basename => $basename, type => 'vcf', metadata => $vcf_meta);
                my $vcf_recal_path = $vcf_recal_file->path;
                
                my $cmd = $self->gatk_prefix($req->memory) . qq[ -T ApplyRecalibration -R $ref --input $vcf_path -recalFile $recal_file_path -tranchesFile $tranches_file_path -o $vcf_recal_path -mode $mode $this_apply_recal_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_apply_recalibration', 'apply_recalibration_and_check', [$cmd, $req, { output_files => [$vcf_recal_file] }]);
                $idx++;
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->determine_gatk_version(),
                    summary => qq[java \$jvm_args -jar GenomeAnalysisTK.jar -T ApplyRecalibration -R \$reference_fasta --input \$vcf_file -recalFile $mode.recal -tranchesFile $mode.recal.tranches -o \$recalibrated_vcf_file -mode $mode $apply_recal_opts]
                )
            );
        };
    }
    
    method outputs_definition {
        return { recalibrated_vcfs => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a recalibrated vcf file for each input vcf') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run GATK ApplyRecalibration to apply VQSLOD scores and recalibration tranche filters to VCF files processed by GATK VariantRecalibrator.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method apply_recalibration_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /-input (\S+) .* -o (\S+) /;
        my $input_file = VRPipe::File->get(path => $input_path);
        
        my $input_records = $input_file->num_records;
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_records = $output_file->num_records;
        
        # Should have extra header lines
        if ($output_records < $input_records) {
            $output_file->unlink;
            $self->throw("Output VCF has $output_records records, less than input $input_records");
        }
        else {
            return 1;
        }
    }
}

1;
