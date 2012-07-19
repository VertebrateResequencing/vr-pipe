
=head1 NAME

VRPipe::Steps::bismark_methylation_extractor - a step

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

class VRPipe::Steps::bismark_methylation_extractor with VRPipe::StepRole {
    use File::Basename;
    use Data::Dumper;
    
    method options_definition {
        return {
            bismark_meth_extractor_exe => VRPipe::StepOption->create(description => 'path to your bismark methylation extractor executable',                          optional => 1, default_value => $ENV{BISMARK_METH_EXR_EXE}),
            paired_end                 => VRPipe::StepOption->create(description => 'Run on files generated with paired end data (default is to run on single end).', optional => 1, default_value => 0)
        
        };
    }
    
    method inputs_definition {
        return { sam_file => VRPipe::StepIODefinition->create(type => 'any', max_files => -1, description => 'A sam file') };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $exe      = $options->{'bismark_meth_extractor_exe'};
            my $paired   = $options->{'paired_end'};
            my ($infile) = @{ $self->inputs->{sam_file} };
            my $name     = $infile->basename;
            
            my $outfile1 = $self->output_file(output_key => 'meth_calls_CHH',
                                              basename   => 'CHH_context_' . $name . '.txt',
                                              type       => 'txt',
                                              metadata   => $infile->metadata);
            
            my $outfile2 = $self->output_file(output_key => 'meth_calls_CpG',
                                              basename   => 'CHG_context_' . $name . '.txt',
                                              type       => 'txt',
                                              metadata   => $infile->metadata);
            
            my $outfile3 = $self->output_file(output_key => 'meth_calls_CHG',
                                              basename   => 'CpG_context_' . $name . '.txt',
                                              type       => 'txt',
                                              metadata   => $infile->metadata);
            my $infile_path = $infile->path;
            
            my $req = $self->new_requirements(memory => 1500, time => 1);
            
            my $cmd = $paired ? $exe . " -p --comprehensive $infile_path" : $exe . " -s --comprehensive $infile_path"; #> $outfile_path";
            $self->dispatch([qq[$cmd], $req, { output_files => [$outfile1, $outfile2, $outfile3] }]);
          }
    }
    
    method outputs_definition {
        return { meth_calls_CHH => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Text file listing methylation calls by position'),
                 meth_calls_CpG => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Text file listing methylation calls by position'),
                 meth_calls_CHG => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Text file listing methylation calls by position') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Step for the Methylation Extractor tool bundled with Bismark to produce a table of methylation calls by position";
    }
    
    method max_simultaneous {
        return 0;                                                                                                      # meaning unlimited
    }
}
