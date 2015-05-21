
=head1 NAME

VRPipe::Steps::bam_to_cram - a step

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

class VRPipe::Steps::bam_to_cram with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe    => VRPipe::StepOption->create(description => 'Path to samtools 1.0 or greater executable',                    optional => 1, default_value => 'samtools'),
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file used to do the mapping', optional => 1),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more coordinate-sorted bam files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $options  = $self->options;
            my $samtools = $options->{samtools_exe};
            my $ref      = $options->{reference_fasta} ? " -T $options->{reference_fasta}" : '';
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools view -C \$bam > \$cram"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_path = $bam->path;
                my $basename = $bam->basename;
                $basename =~ s/bam$/cram/;
                my $cram_file = $self->output_file(
                    output_key => 'cram_files',
                    basename   => $basename,
                    type       => 'cram',
                    metadata   => $bam->metadata
                );
                my $cmd = qq[$samtools view$ref -C $bam_path > ] . $cram_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_to_cram', 'cram_and_check', [$cmd, $req, { output_files => [$cram_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { cram_files => VRPipe::StepIODefinition->create(type => 'cram', max_files => -1, description => 'a cram file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts BAM files to CRAM files using samtools";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method cram_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /-C (\S+) > (\S+)$/;
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
            $out_file->add_metadata({ reads => $actual_reads });
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output CRAM file, yet there were $expected_reads reads in the original BAM file");
        }
    }
}

1;
