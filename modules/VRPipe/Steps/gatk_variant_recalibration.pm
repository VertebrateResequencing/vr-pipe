
=head1 NAME

VRPipe::Steps::gatk_variant_recalibration - a step

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

# java -Xmx4g -jar GenomeAnalysisTK.jar \
#   -T VariantRecalibrator \
#   -R reference/human_g1k_v37.fasta \
#   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.hg19.subset.vcf \
#   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
#   -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
#   -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 dbsnp_132.b37.vcf \
#   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
#   -recalFile path/to/output.recal \
#   -tranchesFile path/to/output.tranches \
#   -rscriptFile path/to/output.plots.R

class VRPipe::Steps::gatk_variant_recalibration extends VRPipe::Steps::gatk {
    around options_definition {
        return {
            %{ $self->$orig },
            variant_recalibration_options => VRPipe::StepOption->create(description => 'options for GATK VariantRecalibrator, including -resource and -annotation options for training the gaussian mixture model, but excluding reference genome, mode, input and output files.'),
            variant_recalibration_mode    => VRPipe::StepOption->create(description => 'mode', optional => 1, default_value => 'BOTH'),
        };
    }
    
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => '1 or more vcf files') };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $recal_opts = $options->{variant_recalibration_options};
            if ($recal_opts =~ /$ref|VariantRecalibrator|-mode|-recalFile|-tranchesFile|-rscriptFile/) {
                $self->throw("variant_recalibration_options should not include the reference, resource file options or VariantRecalibrator task command");
            }
            my $mode = $options->{variant_recalibration_mode};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->determine_gatk_version(),
                    summary => qq[java \$jvm_args -jar GenomeAnalysisTK.jar -T VariantRecalibrator -R \$reference_fasta -I \$vcf_file(s) -recalFile $mode.recal -tranchesFile $mode.recal.tranches -rscriptFile $mode.recal.r -mode $mode $recal_opts]
                )
            );
            
            my $req = $self->new_requirements(memory => 4000, time => 1);
            my $jvm_args = $self->jvm_args($req->memory);
            
            my @vcfs = map { '-input ' . $_->path } @{ $self->inputs->{vcf_files} };
            my $vcf_list = join ' ', @vcfs;
            
            my $recal_file         = $self->output_file(output_key => 'recalibration_file',      basename => "$mode.recal",              type => 'txt');
            my $tranches_file      = $self->output_file(output_key => 'tranches_file',           basename => "$mode.recal.tranches",     type => 'txt');
            my $tranches_plot_file = $self->output_file(output_key => 'tranches_plot_file',      basename => "$mode.recal.tranches.pdf", type => 'any');
            my $rscript_file       = $self->output_file(output_key => 'rscript_file',            basename => "$mode.recal.r",            type => 'txt');
            my $recal_plot_file    = $self->output_file(output_key => 'recalibration_plot_file', basename => "$mode.recal.r.pdf",        type => 'any');
            
            my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T VariantRecalibrator -R $ref ] . $vcf_list . qq[ -recalFile ] . $recal_file->path . qq[ -tranchesFile ] . $tranches_file->path . qq[ -rscriptFile ] . $rscript_file->path . qq[ -mode $mode $recal_opts];
            $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_variant_recalibration', 'variant_recal_and_check', [$this_cmd, $req, { output_files => [$recal_file, $tranches_file, $tranches_plot_file, $rscript_file, $recal_plot_file] }]);
        };
    }
    
    method outputs_definition {
        return {
            recalibration_file      => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'Variant recalibration data file'),
            tranches_file           => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'Variant recalibration tranches file'),
            tranches_plot_file      => VRPipe::StepIODefinition->create(type => 'any', max_files => 1, description => 'Plots to visualise the recalibration data', check_existence => 0),
            rscript_file            => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'Rscript file generated by the VQSR to aid in visualisation of the input data and learned model.'),
            recalibration_plot_file => VRPipe::StepIODefinition->create(type => 'any', max_files => 1, description => 'PDF file associated with rscript file.', check_existence => 0)
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Recalibrates variant calls using Variant Quality Score Recalibration (VQSR)";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method variant_recal_and_check (ClassName|Object $self: Str $cmd_line) {
        # NEED TO COME UP WITH A CHECK HERE
        # my (@input_paths) = $cmd_line =~ /-input (\S+)/g;
        # my $output_path = $cmd_line =~ /-o (\S+)/;
        #
        # my $input_lines = $input_file->lines;
        #
        # $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        #
        # my $output_file = VRPipe::File->get(path => $output_path);
        # $output_file->update_stats_from_disc;
        # my $output_lines = $output_file->lines;
        #
        # # Should have extra header lines
        # unless ($output_lines >= $input_lines) {
        #     $output_file->unlink;
        #     $self->throw("Output VCF has $output_lines lines, less than input $input_lines");
        # }
        # else {
        #     return 1;
        # }
    }
}

1;
