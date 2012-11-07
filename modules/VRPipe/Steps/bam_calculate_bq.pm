
=head1 NAME

VRPipe::Steps::bam_calculate_bq - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bam_calculate_bq with VRPipe::StepRole {
    method options_definition {
        return {
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file used to do the mapping'),
            samtools_exe    => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            ),
            samtools_calmd_options => VRPipe::StepOption->create(
                description   => 'options to samtools calmd',
                optional      => 1,
                default_value => '-Erb'
            )
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more bam files',
            )
        };
    
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $ref = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $samtools   = $options->{samtools_exe};
            my $calmd_opts = $options->{samtools_calmd_options};
            if ($calmd_opts =~ /$ref|calmd/) {
                $self->throw("samtools_calmd_options should not include the reference fasta or the calmd subcommand");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools calmd $calmd_opts \$bam_file \$reference_fasta > \$bq_bam_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 3000, time => 2);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_base = $bam->basename;
                my $bq_base  = $bam_base;
                $bq_base =~ s/bam$/calmd.bam/;
                my $bam_meta = $bam->metadata;
                my $bq_bam_file = $self->output_file(output_key => 'bq_bam_files', basename => $bq_base, type => 'bam', metadata => $bam_meta);
                
                my $this_cmd = "$samtools calmd $calmd_opts " . $bam->path . " $ref > " . $bq_bam_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_calculate_bq', 'calmd_and_check', [$this_cmd, $req, { output_files => [$bq_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bq_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with BQ tag and good NM & MD tags',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Corrects NM & MD tags and calculates BQ (BAQ scores) to aid downstream variant calling";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method calmd_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /(\S+) \S+ > (\S+)/;
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
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output bam file, yet there were $expected_reads reads in the original bam file");
        }
    }
}

1;
