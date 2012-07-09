=head1 NAME

VRPipe::Steps::trimmomatic - a step

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

class VRPipe::Steps::trimmomatic extends VRPipe::Steps::java {
use File::Basename; 
       
    around options_definition {
        return { %{$self->$orig},
                 trimmomatic_jar_path => VRPipe::StepOption->get(description => 'path to Trimmomatic jar file', optional => 1, default_value => "$ENV{TRIMMOMATIC}"),
                 paired_end => VRPipe::StepOption->get(description => 'Run in Paired End mode (default is for single end).', optional => 1, default_value => "0"),
                 trimmomatic_step_options => VRPipe::StepOption->get(description => 'String of the step options for Trimmomatic.', optional => 1, default_value => 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'),
                 log_file => VRPipe::StepOption->get(description => 'Path for log file.', optional => 1, default_value => '')
        };
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->get(type => 'fq', max_files => -1, description => '1 or more fastq files to trim reads.')
        };
    }
    method body_sub {
        return sub { 
           my $self = shift;
           my $options = $self->options;
           $self->handle_standard_options($options);

           my $trimmomatic_jar_path = $options->{'trimmomatic_jar_path'};
           my $end_type = $options->{'paired_end'} ? 'org.usadellab.trimmomatic.TrimmomaticPE' : 'org.usadellab.trimmomatic.TrimmomaticSE';
           my $qual_enc = '-phred33';
           my $req = $self->new_requirements(memory => 1500, time => 1);
           my $jvm_args = $self->jvm_args($req->memory);
                    
           my $log_file = $self->output_file( output_key => 'trimmomatic_files',
                                          basename => 'trimmomatic.log',
                                          type => 'txt',
                                          #metadata => $seq_file->metadata
           )->path;
           
           my $step_options = $options->{'trimmomatic_step_options'};
         
           # IF SINGLE END
           # warn Dumper $self->inputs->{fastq_files}; 
           foreach my $seq_file  (@{$self->inputs->{fastq_files}}) {
               my ($name) = fileparse( $seq_file->basename, ('.fastq') );
               my $out_file = $self->output_file( output_key => 'trimmomatic_files',
                                          basename => $name .'.trim.fastq',
                                          type => 'fq',
                                          metadata => $seq_file->metadata);
               my $out_file_path = $out_file->path; 
               my $seq_file_path = $seq_file->path;
               my $cmd = $self->java_exe . " $jvm_args -classpath $trimmomatic_jar_path org.usadellab.trimmomatic.TrimmomaticSE $qual_enc -trimlog $log_file $seq_file_path $out_file_path $step_options";
              $self->dispatch([ qq[$cmd], $req, { output_files => [ $out_file ] } ]);
           }
           # IF PAIRED END
           # foreach pair of inputs 
           #   $self->dispatch();
       };
    }

    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Step for the Trimmomatic Read Trimmer";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}
