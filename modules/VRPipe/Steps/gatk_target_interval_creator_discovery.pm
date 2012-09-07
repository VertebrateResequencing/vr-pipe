
=head1 NAME

VRPipe::Steps::gatk_target_interval_creator_discovery - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::gatk_target_interval_creator_discovery extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{ $self->$orig }, gatk_target_intervals_discovery_options => VRPipe::StepOption->create(description => 'command line options for GATK RealignerTargetCreator, when used to discover locations to realign around, excluding -R and -T', optional => 1) };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more bam files'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $intervals_opts = $options->{gatk_target_intervals_discovery_options};
            if ($intervals_opts =~ /$ref|RealignerTargetCreator/) {
                $self->throw("gatk_target_intervals_discovery_options should not include the reference or RealignerTargetCreator task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference_fasta -I $input_bam -o $intervals_file ' . $intervals_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 1);
            my $jvm_args = $self->jvm_args($req->memory);
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_base       = $bam->basename;
                my $intervals_file = $self->output_file(
                    output_key => 'intervals_file',
                    basename   => qq[$bam_base.intervals],
                    type       => 'txt',
                    metadata   => { source_bam => $bam->path->stringify }
                );
                
                my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T RealignerTargetCreator -R $ref -o ] . $intervals_file->path . ' -I ' . $bam->path . ' ' . $intervals_opts;
                $self->dispatch([$this_cmd, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            intervals_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'GATK intervals file for known indel sites',
                metadata    => { source_bam => 'the bam that was used to discover locations to realign around' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Creates target intervals file for a bam, discovering indels that could be realigned around.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
