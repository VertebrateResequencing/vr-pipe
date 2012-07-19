
=head1 NAME

VRPipe::Steps::bam_index - a step

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

class VRPipe::Steps::bam_index with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->create(description   => 'path to your samtools executable',
                                                            optional      => 1,
                                                            default_value => 'samtools') };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to index') };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $samtools = $options->{samtools_exe};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_path = $bam->path;
                my $bai_file = $self->output_file(output_key => 'bai_files',
                                                  output_dir => $bam->dir,
                                                  basename   => $bam->basename . '.bai',
                                                  type       => 'bin',
                                                  metadata   => $bam->metadata);
                my $bai_path = $bai_file->path;
                my $cmd      = qq[$samtools index $bam_path $bai_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_index', 'index_and_check', [$cmd, $req, { output_files => [$bai_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { bai_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'a .bai file for each input bam file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Indexes bam files using samtools";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method index_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $bai_path) = $cmd_line =~ /index (\S+) (\S+)/;
        $bam_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $bai_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $bam_file = VRPipe::File->get(path => $bam_path);
        my $bai_file = VRPipe::File->get(path => $bai_path);
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $bai_file->update_stats_from_disc(retries => 3);
        
        my $correct_magic = [qw(102 101 111 001)];
        
        if ($bai_file->check_magic($bai_file->path, $correct_magic)) {
            return 1;
        }
        else {
            $bai_file->unlink;
            $self->throw("cmd [$cmd_line] failed because index file $bai_path had the incorrect magic");
        }
    }
}

1;
