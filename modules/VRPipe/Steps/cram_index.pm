
=head1 NAME

VRPipe::Steps::cram_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012-2015 Genome Research Limited.

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

class VRPipe::Steps::cram_index with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->create(description => 'Path to samtools 1.0 or greater executable', optional => 1, default_value => 'samtools'), };
    }
    
    method inputs_definition {
        return { cram_files => VRPipe::StepIODefinition->create(type => 'cram', max_files => -1, description => '1 or more cram files') };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $samtools = $options->{samtools_exe};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools index \$cram"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $cram (@{ $self->inputs->{cram_files} }) {
                my $cram_path = $cram->path;
                my $crai_file = $self->output_file(
                    output_key => 'cram_index_files',
                    output_dir => $cram->dir,
                    basename   => $cram->basename . '.crai',
                    type       => 'bin',
                    metadata   => $cram->metadata
                );
                my $cmd = qq[$samtools index $cram_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::cram_index', 'cram_index_and_check', [$cmd, $req, { output_files => [$crai_file], block_and_skip_if_ok => 1 }]);
            }
        };
    }
    
    method outputs_definition {
        return { cram_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'a cram index file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Indexes cram files using samtools";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method cram_index_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($cram_path) = $cmd_line =~ /index (\S+)/;
        $cram_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        my $crai_path = $cram_path . '.crai';
        
        my $cram_file = VRPipe::File->get(path => $cram_path);
        my $crai_file = VRPipe::File->get(path => $crai_path);
        
        $crai_file->update_stats_from_disc(retries => 3);
        
        my $crai_already_there = 0;
        if ($cram_file->e) {
            my $cram_timestamp = (stat $cram_path)[9];
            my $crai_timestamp = (stat $crai_path)[9];
            if ($crai_timestamp >= $cram_timestamp) {
                $crai_already_there = 1;
            }
            else {
                $crai_file->remove;
            }
        }
        
        unless ($crai_already_there) {
            $cram_file->disconnect;
            system($cmd_line) && $self->throw("failed to run [$cmd_line]");
            $crai_file->update_stats_from_disc(retries => 3);
        }
        my $correct_magic = [qw(037 213 010 000)];
        
        if ($crai_file->check_magic($crai_file->path, $correct_magic)) {
            return 1;
        }
        else {
            $crai_file->unlink;
            $self->throw("cmd [$cmd_line] failed because index file $crai_path had the incorrect magic");
        }
    }
}

1;
