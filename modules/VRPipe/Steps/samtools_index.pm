
=head1 NAME

VRPipe::Steps::samtools_index - a step

=head1 DESCRIPTION

*** more documentation to come

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

class VRPipe::Steps::samtools_index with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            )
        };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'aln', max_files => -1, description => '1 or more BAM or CRAM files to index') };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $samtools = $options->{samtools_exe};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $aln (@{ $self->inputs->{bam_files} }) {
                my $aln_path   = $aln->path;
                my $suffix     = $aln->type eq 'bam' ? '.bai' : '.crai';
                my $index_file = $self->output_file(
                    output_key => 'index_files',
                    output_dir => $aln->dir,
                    basename   => $aln->basename . $suffix,
                    type       => 'bin',
                    metadata   => $aln->metadata
                );
                my $index_path = $index_file->path;
                my $cmd        = qq[$samtools index $aln_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::samtools_index', 'index_and_check', [$cmd, $req, { output_files => [$index_file], block_and_skip_if_ok => 1 }]);
            }
        };
    }
    
    method outputs_definition {
        return { index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'a .bai or .crai file for each input BAM/CRAM file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Indexes BAM and CRAM files using samtools";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method index_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($aln_path) = $cmd_line =~ /index (\S+)/;
        $aln_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        my $aln_file = VRPipe::File->get(path => $aln_path);
        
        my $suffix = $aln_file->type eq 'bam' ? '.bai' : '.crai';
        my $idx_path = $aln_path . $suffix;
        
        my $idx_file = VRPipe::File->get(path => $idx_path);
        
        $idx_file->update_stats_from_disc(retries => 3);
        
        my $idx_already_there = 0;
        if ($idx_file->e) {
            my $aln_timestamp = (stat $aln_path)[9];
            my $idx_timestamp = (stat $idx_path)[9];
            if ($idx_timestamp >= $aln_timestamp) {
                $idx_already_there = 1;
            }
            else {
                $idx_file->remove;
            }
        }
        
        unless ($idx_already_there) {
            $aln_file->disconnect;
            system($cmd_line) && $self->throw("failed to run [$cmd_line]");
            $idx_file->update_stats_from_disc(retries => 3);
        }
        
        if ($suffix eq 'bai') {
            my $bai = VRPipe::FileType->create('hts', { file => $idx_path });
            if ($bai->hts_file_type =~ /^BAI/) {
                return 1;
            }
            else {
                $idx_file->unlink;
                $self->throw("cmd [$cmd_line] failed because index file $idx_path had the incorrect magic");
            }
        }
        else {
            return 1;
        }
    }
}

1;
