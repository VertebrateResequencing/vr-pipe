
=head1 NAME

VRPipe::Steps::bam_shuffle_by_name - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::bam_shuffle_by_name with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            ),
            bamshuf_options => VRPipe::StepOption->create(description => 'command line options for samtools bamshuf, excluding the -O option', optional => 1)
        };
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
            
            my $samtools = $options->{samtools_exe};
            my $opts     = $options->{bamshuf_options};
            if ($opts =~ /-O/) {
                $self->throw("bamshuf_options should not include the -O option");
            }
            
            my $req = $self->new_requirements(memory => 1000, time => 2);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $in_base  = $bam->basename;
                my $out_base = $in_base;
                $out_base =~ s/\.bam$/.shuf.bam/;
                my $shuf_bam_file = $self->output_file(output_key => 'name_shuffled_bam_files', basename => $out_base, type => 'bam', metadata => $bam->metadata);
                
                my $out_prefix = $shuf_bam_file->path;
                $out_prefix =~ s/\.bam$//;
                my $this_cmd = "$samtools bamshuf $opts " . $bam->path . " $out_prefix";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_shuffle_by_name', 'shuffle_and_check', [$this_cmd, $req, { output_files => [$shuf_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            name_shuffled_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with alignments grouped by name'
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Shuffle and group alignments by name";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method shuffle_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /(\S+) (\S+)$/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path .= '.bam';
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->_filetype->check_records_vs_input($in_file, $cmd_line);
    }
}

1;
