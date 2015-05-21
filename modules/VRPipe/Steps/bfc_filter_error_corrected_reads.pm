
=head1 NAME

VRPipe::Steps::bfc_filter_error_corrected_reads - a step

=head1 DESCRIPTION

Apply filtering to bfc error corrected reads

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

class VRPipe::Steps::bfc_filter_error_corrected_reads with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return {
            bfc_exe            => VRPipe::StepOption->create(description => 'path to bfc executable',                  default_value => 'bfc'),
            bfc_filter_options => VRPipe::StepOption->create(description => 'options for bcf error correction filter', default_value => '-1s 3g -k 61 -t 16'),
        };
    }
    
    method inputs_definition {
        return {
            bfc_error_corrected_reads => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'bfc error corrected fastq reads',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $bfc      = $options->{bfc_exe};
            my $bfc_opts = $options->{bfc_filter_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bfc',
                    version => VRPipe::StepCmdSummary->determine_version(qq[$bfc -v], '^(.+)$'),
                    summary => qq[bfc $bfc_opts \$output.ec.fq.gz | gzip -1 > \$output.flt.fq.gz]
                )
            );
            
            my ($cpus) = $bfc_opts =~ m/-t\s*(\d+)/;
            my $req = $self->new_requirements(memory => 40000, time => 1, $cpus ? (cpus => $cpus) : ());
            foreach my $fq (@{ $self->inputs->{bfc_error_corrected_reads} }) {
                my $filt_fq = $self->output_file(
                    output_key => 'bfc_filtered_error_corrected_reads',
                    basename   => $fq->basename,
                    type       => 'fq',
                    metadata   => $fq->metadata
                );
                my $this_cmd = qq[$bfc $bfc_opts ] . $fq->path . qq[ | gzip -1 > ] . $filt_fq->path;
                $self->dispatch([$this_cmd, $req, { output_files => [$filt_fq] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bfc_filtered_error_corrected_reads => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'bfc filtered error corrected fastq reads',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Apply filtering to bfc error corrected reads";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
