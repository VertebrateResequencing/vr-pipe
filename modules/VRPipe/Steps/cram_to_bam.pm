
=head1 NAME

VRPipe::Steps::cram_to_bam - a step

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

class VRPipe::Steps::cram_to_bam with VRPipe::StepRole {
    method options_definition {
        return { samtools_exe => VRPipe::StepOption->create(description => 'Path to samtools 1.0 or greater executable', optional => 1, default_value => 'samtools'), };
    }
    
    method inputs_definition {
        return {
            cram_files => VRPipe::StepIODefinition->create(type => 'cram', max_files => -1, description => '1 or more cram files'),
        };
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
                    summary => "samtools view -b \$cram > \$bam"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $cram (@{ $self->inputs->{cram_files} }) {
                my $cram_path = $cram->path;
                my $basename  = $cram->basename;
                $basename =~ s/cram$/bam/;
                my $bam_file = $self->output_file(
                    output_key => 'bam_files',
                    basename   => $basename,
                    type       => 'bam',
                    metadata   => $cram->metadata,
                );
                my $cmd = qq[$samtools view -b $cram_path > ] . $bam_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::cram_to_bam', 'cram_to_bam_and_check', [$cmd, $req, { output_files => [$bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => 'a bam file'), };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts CRAM files to BAM files using samtools";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method cram_to_bam_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /-b (\S+) > (\S+)$/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output BAM file, yet there were $expected_reads reads in the original CRAM file");
        }
    }
}

1;
