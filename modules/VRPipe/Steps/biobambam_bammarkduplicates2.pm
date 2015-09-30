
=head1 NAME

VRPipe::Steps::biobambam_bammarkduplicates2 - a step

=head1 DESCRIPTION

Runs biobambam's bammarkduplicates2 to mark duplicates.

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

class VRPipe::Steps::biobambam_bammarkduplicates2 with VRPipe::StepRole {
    method options_definition {
        return {
            bammarkduplicates2_exe     => VRPipe::StepOption->create(description => 'path to bammarkduplicates2 executable',                                             optional => 1, default_value => 'bammarkduplicates2'),
            bammarkduplicates2_options => VRPipe::StepOption->create(description => 'bammarkduplicates2 options (excluding arguments that set input/output file names)', optional => 1, default_value => 'resetdupflag=1 outputthreads=4'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',                                          # cram or bam
                max_files   => -1,
                description => '1 or more coordinate sorted BAM or CRAM files',
                metadata    => {
                    bases    => 'total number of base pairs',
                    reads    => 'total number of reads (sequences)',
                    paired   => '0=single ended reads only; 1=paired end reads present',
                    optional => ['bases', 'reads', 'paired']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                    = shift;
            my $options                 = $self->options;
            my $bammarkduplicates2_exe  = $options->{bammarkduplicates2_exe};
            my $bammarkduplicates2_opts = $options->{bammarkduplicates2_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bammarkduplicates2',
                    version => VRPipe::StepCmdSummary->determine_version($bammarkduplicates2_exe . ' --version', '^This is biobambam\d? version (.+)\.$'),
                    summary => "bammarkduplicates2 $bammarkduplicates2_opts O=\$markdup_file index=1 indexfilename=\$markdup_index M=\$metrics_file"
                )
            );
            
            my ($input_threads)  = $bammarkduplicates2_opts =~ m/inputthreads=(\d+)/;
            my ($output_threads) = $bammarkduplicates2_opts =~ m/outputthreads=(\d+)/;
            my $cpus             = 1;
            $cpus = $input_threads  if ($input_threads  && $input_threads > $cpus);
            $cpus = $output_threads if ($output_threads && $output_threads > $cpus);
            my $req = $self->new_requirements(memory => 8000, time => 1, cpus => $cpus);
            
            foreach my $aln (@{ $self->inputs->{bam_files} }) {
                my $basename = $aln->basename;
                $basename =~ s/\.(cr|b)am$//;
                my $markdup_file = $self->output_file(
                    output_key => 'markdup_files',
                    basename   => $basename . '.bam',
                    type       => 'bam',
                    metadata   => $aln->metadata
                );
                my $markdup_index_file = $self->output_file(
                    output_key => 'markdup_index_files',
                    basename   => $basename . '.bam.bai',
                    type       => 'bai',
                    metadata   => $aln->metadata
                );
                my $markdup_metrics_file = $self->output_file(
                    output_key => 'markdup_metrics_files',
                    basename   => $basename . '.metrics',
                    type       => 'txt',
                    metadata   => $aln->metadata
                );
                my $markdup_path = $markdup_file->path;
                my $this_cmd     = "$bammarkduplicates2_exe $bammarkduplicates2_opts I=" . $aln->path . " O=$markdup_path index=1 indexfilename=" . $markdup_index_file->path . " M=" . $markdup_metrics_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::biobambam_bammarkduplicates2', 'markdup_and_check', [$this_cmd, $req, { output_files => [$markdup_file, $markdup_index_file, $markdup_metrics_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            markdup_files         => VRPipe::StepIODefinition->create(type => 'aln', max_files => -1, description => 'a BAM file with duplicates marked'),
            markdup_index_files   => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'BAI index file for merged, markdup-ed BAM file'),
            markdup_metrics_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'a text file with metrics from biobambam duplicate marking'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Marks duplicates for BAM or CRAM files using biobambam bammarkduplicates2";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method markdup_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /I=(\S+).+O=(\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->_filetype->check_records_vs_input($in_file, $cmd_line);
    }

}

1;
