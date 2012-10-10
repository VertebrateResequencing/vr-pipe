
=head1 NAME

VRPipe::Steps::gatk_variant_recalibration_for_snps - a step

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

class VRPipe::Steps::gatk_variant_recalibration_for_snps extends VRPipe::Steps::gatk_variant_recalibration {
    around options_definition {
        my $opts_def = $self->$orig;
        delete $opts_def->{variant_recalibration_options};
        delete $opts_def->{variant_recalibration_mode};
        return { %{$opts_def}, snp_recalibration_options => VRPipe::StepOption->create(description => 'Command line options for GATK VariantRecalibrator when in SNP mode') };
    }
    
    around options {
        my $options = $self->$orig();
        $options->{variant_recalibration_options} = $options->{snp_recalibration_options};
        $options->{variant_recalibration_mode}    = 'SNP';
        return $options;
    
    }
    
    method description {
        return "Recalibrates SNP calls using Variant Quality Score Recalibration (VQSR)";
    }
}

1;
