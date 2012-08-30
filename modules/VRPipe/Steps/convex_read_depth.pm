
=head1 NAME

VRPipe::Steps::convex_read_depth - a step

=head1 DESCRIPTION

Runs the java ReadDepth class from the CoNVex package, generating region read
depths from a bam

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::convex_read_depth extends VRPipe::Steps::java {
    around options_definition {
        return {
            %{ $self->$orig },
            'convex_classpath' => VRPipe::StepOption->create(description => 'convex classpath to all convex package jars'),
            'regions_file'     => VRPipe::StepOption->create(description => 'regions file for which to generate read depths'),
            'chr_prefix'       => VRPipe::StepOption->create(description => 'chromosome name prefix within the bam', optional => 1),
        };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to call variants') };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $convex_classpath = $options->{'convex_classpath'};
            my $regions_file     = $options->{'regions_file'};
            my $chr_prefix       = $options->{'chr_prefix'};
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            my $jvm_args = $self->jvm_args($req->memory);
            
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $bam_path = $bam_file->path;
                my $basename = $bam_file->basename;
                $basename =~ s/\.bam$/.rd.txt/;
                
                my $rd_file = $self->output_file(
                    output_key => 'rd_files',
                    basename   => $basename,
                    type       => 'txt',
                    metadata   => { source_bam => $bam_file->path->stringify, batch => 1 }
                ); # batch metadata is used to merge up for L2R pipeline
                my $rd_path = $rd_file->path;
                
                my $cmd = $self->java_exe . " $jvm_args -classpath $convex_classpath ReadDepth -bam_file $bam_path -regions_file $regions_file";
                $cmd .= " -chr_prefix $chr_prefix" if $chr_prefix;
                $cmd .= "  -rd_file $rd_path;";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::convex_read_depth', 'run_read_depth', [$cmd, $req, { output_files => [$rd_file] }]);
            
            }
        };
    
    }
    
    method outputs_definition {
        return { rd_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'a read depths file for each input bam'), };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs CoNVex ReadDepth, generating region read depths text file from a bam";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method run_read_depth (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ /-rd_file (\S+);$/;
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        unless ($output_lines > 0) {
            $output_file->unlink;
            $self->throw("Output read depth file is zero length");
        }
        else {
            return 1;
        }
    
    }
}

1;
