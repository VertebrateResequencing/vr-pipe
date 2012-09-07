
=head1 NAME

VRPipe::Steps::bam_name_sort - a step

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

class VRPipe::Steps::bam_name_sort with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            ),
            samtools_sort_options => VRPipe::StepOption->create(description => 'command line options for samtools sort, excluding -n and -o', optional => 1)
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
            
            my $opts      = "sort -n";
            my $user_opts = $options->{samtools_sort_options};
            my $mem       = 768;                              # the samtools sort default
            if ($user_opts) {
                $user_opts =~ s/-o//;
                $opts .= ' ' . $user_opts;
                if ($user_opts =~ /-m (\d+)([KMG])/) {
                    $mem = $1;
                    my $m_unit = $2;
                    if ($m_unit eq 'K') {
                        $mem /= 1000;
                    }
                    elsif ($m_unit eq 'G') {
                        $mem *= 1000;
                    }
                }
            }
            
            my $req = $self->new_requirements(memory => $mem == 768 ? 3000 : (4 * $mem), time => 2); # in some cases samtools can use way more than the -m specified, and is very segfault happy
            my $memory = $req->memory;
            
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $in_base  = $bam->basename;
                my $out_base = $in_base;
                $out_base =~ s/\.bam$//;
                $out_base .= '.name_sorted.bam';
                my $bam_meta = $bam->metadata;
                my $sort_bam_file = $self->output_file(output_key => 'name_sorted_bam_files', basename => $out_base, type => 'bam', metadata => $bam_meta);
                
                my $out_prefix = $sort_bam_file->path;
                $out_prefix =~ s/\.bam$//;
                my $this_cmd = "$samtools $opts " . $bam->path . " $out_prefix";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_name_sort', 'sort_and_check', [$this_cmd, $req, { output_files => [$sort_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            name_sorted_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a name-sorted bam file'
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Name-sorts a bam file, producing a new bam file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method sort_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /(\S+) (\S+)$/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path .= '.bam';
        
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
