
=head1 NAME

VRPipe::Steps::gatk_base_recalibrator - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

# java -Xmx4g -jar GenomeAnalysisTK.jar \
#   -T BaseRecalibrator \
#   -I my_reads.bam \
#   -R resources/Homo_sapiens_assembly18.fasta \
#   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
#   -knownSites another/optional/setOfSitesToMask.vcf \
#   -o recal_data.grp

class VRPipe::Steps::gatk_base_recalibrator extends VRPipe::Steps::gatk_v2 {
    around options_definition {
        return {
            %{ $self->$orig },
            base_recalibrator_options => VRPipe::StepOption->create(description => 'command line options for GATK BaseRecalibrator -- must include -cov options; excludes the -knownSites option(s) which are set by another StepOption', optional => 1, default_value => '-l INFO -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y -L MT -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate'),
        };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files') };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = $options->{reference_fasta};
            
            my $recal_options = $options->{base_recalibrator_options};
            if ($recal_options =~ /$ref|-I |--input_file|-o | --output|BaseRecalibrator/) {
                $self->throw("base_recalibrator_options should not include the reference, input files, output files option or BaseRecalibrator task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $reference_fasta -I $bam_file -o $bam_file.recal_data.grp ' . $recal_options
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $recal_base = $bam->basename;
                $recal_base =~ s/bam$/recal_data.grp/;
                my $recal_file = $self->output_file(
                    output_key => 'bam_recalibration_files',
                    basename   => $recal_base,
                    type       => 'txt',
                    metadata   => { source_bam => $bam->path->stringify }
                );
                
                my $temp_dir = $options->{tmp_dir} || $recal_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T BaseRecalibrator -R $ref -I ] . $bam->path . qq[ -o ] . $recal_file->path . qq[ $recal_options];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_base_recalibrator', 'bqsr_and_check', [$this_cmd, $req, { output_files => [$recal_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bam_recalibration_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'recalibration file from GATK BaseRecalibrator',
                metadata    => { source_bam => 'the bam file used to generate this recal file' }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run Base Quality Score Recalibration (BQSR) on a bam file creating a file to be used to recalibrate quality scores";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method bqsr_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($recal_path) = $cmd_line =~ /-o (\S+)/;
        $recal_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $recal_file = VRPipe::File->get(path => $recal_path);
        
        $recal_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $recal_file->update_stats_from_disc(retries => 3);
        
        ## NEED A WAY TO TEST THIS FILE OKAY? CountCovariate use to have and EOF
        # my $eof = `tail -1 $recal_path`;
        # chomp $eof;
        
        # if ($eof eq "EOF") {
        return 1;
        # }
        # else {
        # $recal_file->unlink;
        # $self->throw("cmd [$cmd_line] there was no EOF line at the end of recalibration file $recal_path");
        # }
    }
}

1;
