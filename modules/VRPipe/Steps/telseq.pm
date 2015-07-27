
=head1 NAME

VRPipe::Steps::telseq - a step

=head1 DESCRIPTION

Runs telseq on BAM files to estimate telomere length.

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

class VRPipe::Steps::telseq with VRPipe::StepRole {
    method options_definition {
        return {
            telseq_exe     => VRPipe::StepOption->create(description => 'path to telseq executable', optional => 1, default_value => 'telseq'),
            telseq_options => VRPipe::StepOption->create(description => 'options for telseq',        optional => 1, default_value => ''),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => '1 or more coordinate sorted BAM files',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $options     = $self->options;
            my $telseq      = $options->{telseq_exe};
            my $telseq_opts = $options->{telseq_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'telseq',
                    version => VRPipe::StepCmdSummary->determine_version($telseq, '^Version: (.+)$'),
                    summary => "telseq $telseq_opts \$bam(s) > \$telseq_stats"
                )
            );
            
            my $inputs = $self->inputs->{bam_files};
            my @input_bam_paths;
            if (@$inputs > 5) {
                my $fofn = $self->output_file(basename => "teleseq_input.list", type => 'txt', temporary => 1);
                $fofn->create_fofn($inputs);
                @input_bam_paths = ('-f ' . $fofn->path);
            }
            else {
                @input_bam_paths = map { $_->path } @$inputs;
            }
            
            my $telseq_stats = $self->output_file(
                output_key => 'telseq_stats_files',
                basename   => 'telseq_stats.txt',
                type       => 'txt',
            );
            my $cmd = qq[$telseq $telseq_opts @input_bam_paths > ] . $telseq_stats->path;
            $self->dispatch([$cmd, $self->new_requirements(memory => 1000, time => 1)]);
        };
    }
    
    method outputs_definition {
        return {
            telseq_stats_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'telomere stats from telseq'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs telseq on BAM files to estimate telomere length";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
