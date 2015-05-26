
=head1 NAME

VRPipe::Steps::bfc_error_correct - a step

=head1 DESCRIPTION

Converts BAM or CRAM files to fastq and error corrects reads using bfc

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

class VRPipe::Steps::bfc_error_correct with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe                 => VRPipe::StepOption->create(description => 'path to samtools executable',                                                                      optional      => 1, default_value => 'samtools'),
            samtools_view_options        => VRPipe::StepOption->create(description => 'options for samtools view excluding -u; default set to filter QCFAIL,SECONDARY,DUP,SUPPLEMENTARY', optional      => 1, default_value => '-F0xf00'),
            samtools_bam2fq_options      => VRPipe::StepOption->create(description => 'options for samtools bam2fq',                                                                      optional      => 1, default_value => ''),
            bfc_exe                      => VRPipe::StepOption->create(description => 'path to bfc executable',                                                                           default_value => 'bfc'),
            bfc_error_correction_options => VRPipe::StepOption->create(description => 'options for bcf error correction',                                                                 default_value => '-s 3g -t 16'),
        };
    }
    
    method inputs_definition {
        return { aln_files => VRPipe::StepIODefinition->create(type => 'aln', max_files => -1, description => 'one or more BAM or CRAM files') };
    }
    
    method body_sub {
        return sub {
            my $self                 = shift;
            my $options              = $self->options;
            my $samtools             = $options->{samtools_exe};
            my $samtools_view_opts   = $options->{samtools_view_options};
            my $samtools_bam2fq_opts = $options->{samtools_bam2fq_options};
            my $bfc                  = $options->{bfc_exe};
            my $bfc_opts             = $options->{bfc_error_correction_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools bam2fq $samtools_bam2fq_opts \$input_file(s)"
                )
            );
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bfc',
                    version => VRPipe::StepCmdSummary->determine_version(qq[$bfc -v], '^(.+)$'),
                    summary => qq[bfc $bfc_opts - | gzip -1 > \$output.ec.fq.gz]
                )
            );
            
            my ($cpus) = $bfc_opts =~ m/-t\s*(\d+)/;
            my $req = $self->new_requirements(memory => 60000, time => 1, $cpus ? (cpus => $cpus) : ());
            foreach my $aln (@{ $self->inputs->{aln_files} }) {
                my $prefix = $aln->basename;
                $prefix =~ s/\.(cr|b)am$//;
                my $ec_fastq_file = $self->output_file(
                    output_key => 'bfc_error_corrected_reads',
                    basename   => $prefix . '.fq.gz',
                    type       => 'fq',
                    metadata   => $aln->metadata
                );
                my $bam2fq_cmd;
                if ($samtools_view_opts) {
                    $bam2fq_cmd = qq[$samtools view -u $samtools_view_opts ] . $aln->path . qq[ | $samtools bam2fq - ];
                }
                else {
                    $bam2fq_cmd = qq[$samtools bam2fq $samtools_bam2fq_opts ] . $aln->path;
                }
                my $cmd = qq[$bfc $bfc_opts <($bam2fq_cmd) <($bam2fq_cmd) | gzip -1 > ] . $ec_fastq_file->path;
                $self->dispatch([$cmd, $req, { output_files => [$ec_fastq_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bfc_error_corrected_reads => VRPipe::StepIODefinition->create(type => 'fq', max_files => 1, description => 'a fastq file containing bfc error corrected reads'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts BAM or CRAM files to fastq and error corrects reads using bfc";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
