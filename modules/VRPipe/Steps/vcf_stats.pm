
=head1 NAME

VRPipe::Steps::vcf_stats - a step

=head1 DESCRIPTION

Generate stats file using vcf-stats for each input VCF

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::vcf_stats with VRPipe::StepRole {
    method options_definition {
        return {
            'vcf-stats_options' => VRPipe::StepOption->create(description => 'vcf-stats options'),
            'vcf-stats_exe'     => VRPipe::StepOption->create(
                description   => 'path to vcf-stats executable',
                optional      => 1,
                default_value => 'vcf-stats'
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'vcf files',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options    = $self->options;
            my $stats_exe  = $options->{'vcf-stats_exe'};
            my $stats_opts = $options->{'vcf-stats_options'};
            my $cat_exe;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                my $basename = $vcf_file->basename;
                my $cat_exe = $basename =~ /\.vcf.gz$/ ? 'zcat' : 'cat';
                $basename .= '.stats';
                
                my $stats_file = $self->output_file(output_key => 'stats_file', basename => $basename, type => 'txt', metadata => $vcf_file->metadata);
                
                my $input_path  = $vcf_file->path;
                my $output_path = $stats_file->path;
                
                my $this_cmd = "$cat_exe $input_path | $stats_exe $stats_opts > $output_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_stats', 'vcf_stats', [$this_cmd, $req, { output_files => [$stats_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            stats_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a vcf stats file',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generate stats file for input VCFs";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method vcf_stats (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /^\S+ (\S+) .* (\S+)$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        unless ($output_lines > 0) {
            $output_file->unlink;
            $self->throw("no data in output stats file");
        }
        else {
            return 1;
        }
    }
}

1;
