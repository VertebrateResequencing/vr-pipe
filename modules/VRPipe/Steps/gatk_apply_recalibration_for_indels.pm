
=head1 NAME

VRPipe::Steps::gatk_apply_recalibration_for_indels - a step

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

class VRPipe::Steps::gatk_apply_recalibration_for_indels extends VRPipe::Steps::gatk_apply_recalibration {
    around options_definition {
        my $opts_def = $self->$orig;
        delete $opts_def->{apply_recalibration_options};
        return { %{$opts_def}, apply_recalibration_options_for_indels => VRPipe::StepOption->create(description => 'command line options for GATK ApplyRecalibration for INDELs') };
    }
    
    around options {
        my $options = $self->$orig();
        $options->{apply_recalibration_options} = $options->{apply_recalibration_options_for_indels};
        $options->{apply_recalibration_mode}    = 'INDEL';
        return $options;
    
    }
    
    method description {
        return "Applies recalibration to INDEL calls using Variant Quality Score Recalibration (VQSR)";
    }
}

1;
