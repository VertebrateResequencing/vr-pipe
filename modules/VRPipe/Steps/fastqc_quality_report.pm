
=head1 NAME

VRPipe::Steps::fastqc_quality_report - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

NJWalker <nw11@sanger.ac.uk>.

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

class VRPipe::Steps::fastqc_quality_report with VRPipe::StepRole {
    use File::Basename;
    use File::Which;
    
    method options_definition {
        return { fastqc_exe => VRPipe::StepOption->create(description => 'path to your fastqc executable', optional => 1, default_value => 'fastqc') }
    
    }
    
    method inputs_definition {
        return {
            # sequence file - fastq for now
            fastq_files => VRPipe::StepIODefinition->create(type => 'fq', max_files => -1, description => '1 or more fastq files to calculate quality reports for')
        
        };
    }
    
    method body_sub {
        return sub {
            my $self                             = shift;
            my $options                          = $self->options;
            my $fastqc                           = $options->{fastqc_exe};
            my $fastqc_exe_path                  = which($fastqc) unless file($fastqc)->is_absolute;
            my $fastqc_exe_path_with_interpreter = "perl $fastqc_exe_path";
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => $fastqc_exe_path_with_interpreter, version => VRPipe::StepCmdSummary->determine_version($fastqc . ' --version', '^FastQC v(.+)$'), summary => 'fastqc --noextract file1 '));
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            #my @a =  @{ $self->inputs->{fastq_files} };
            #warn "here @a";
            
            foreach my $seq_file (@{ $self->inputs->{fastq_files} }) {
                my ($name) = fileparse($seq_file->basename, ('.fastq', '.fastq.gz'));
                my $report_file = $self->output_file(
                    output_key => 'fastq_report_file',
                    basename   => $name . '_fastqc.zip',
                    type       => 'bin',
                    metadata   => $seq_file->metadata
                );
                my $seq_file_path   = $seq_file->path;
                my $report_file_dir = $report_file->path->dir;
                my $cmd             = qq[$fastqc_exe_path_with_interpreter --noextract $seq_file_path --outdir $report_file_dir];
                $self->dispatch([$cmd, $req, { output_files => [$report_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { fastq_report_file => VRPipe::StepIODefinition->create(type => 'bin', description => 'a zip file containing the fastqc quality report files') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Produces quality report using fastqc";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}
