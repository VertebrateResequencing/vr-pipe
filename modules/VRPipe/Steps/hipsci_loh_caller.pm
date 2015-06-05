
=head1 NAME

VRPipe::Steps::hipsci_loh_caller - a step

=head1 DESCRIPTION

This step runs hipsci_loh_caller.pl, which calls loss of heterozygosity on
genotype array data in VCF format.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::hipsci_loh_caller with VRPipe::StepRole {
    use VRPipe::Schema;
    
    method options_definition {
        return {
            hipsci_loh_caller_exe => VRPipe::StepOption->create(
                description   => 'path to your hipsci_loh_caller.pl script',
                optional      => 1,
                default_value => 'hipsci_loh_caller.pl'
            ),
            hipsci_loh_caller_options => VRPipe::StepOption->create(
                description   => 'options to hipsci_loh_caller.pl, excluding -l, -v and -c',
                optional      => 1,
                default_value => '-sb'
            ),
            control_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key to check on the input files to extract the sample identifier of the control from',
                default_value => 'sample_control'
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => '1 or more VCF files'
            )
        };
    
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $cmk         = $options->{control_metadata_key};
            my $caller      = $options->{hipsci_loh_caller_exe};
            my $caller_opts = $options->{hipsci_loh_caller_options};
            if ($caller_opts =~ /\s-[lvc]\s/) {
                $self->throw("hipsci_loh_caller_options should not include -l, -v or -c");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'hipsci_loh_caller.pl',
                    version => 0,
                    summary => "hipsci_loh_caller.pl -v \$vcf_file -c \$sample $caller_opts > \$output_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $base = $vcf->basename;
                $base =~ s/\.gz$//;
                $base =~ s/\.vcf$//;
                $base .= '.txt';
                my $meta = $vcf->metadata;
                
                my $sample = $meta->{$cmk};
                $self->throw($vcf->path . " lacks a sample value for the $cmk metadata key") unless $sample;
                
                my $output_file = $self->output_file(output_key => 'loh_caller_output_files', basename => $base,          type => 'txt', metadata  => $meta);
                my $log_file    = $self->output_file(output_key => 'loh_caller_temp_files',   basename => $base . '.log', type => 'txt', temporary => 1);
                
                my $this_cmd = "$caller -v " . $vcf->path . " -c $sample -l " . $log_file->path . " $caller_opts > " . $output_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::hipsci_loh_caller', 'call_and_check', [$this_cmd, $req, { output_files => [$output_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            loh_caller_output_files => VRPipe::StepIODefinition->create(
                type            => 'txt',
                max_files       => -1,
                description     => 'output from the hipsci loh caller',
                check_existence => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Make loss of heterozyosity calls on a VCF using hipsci_loh_caller.pl.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method call_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $log_path, $out_path) = $cmd_line =~ /-v (\S+) .+? -l (\S+) .* > (\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $log_file = VRPipe::File->get(path => $log_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        # check the out file has as many records in as the log files says there
        # should be
        $log_file->update_stats_from_disc(retries => 3);
        $out_file->update_stats_from_disc(retries => 3);
        
        my $lfh = $log_file->openr;
        my $expected_records;
        while (<$lfh>) {
            if (/^(\d+) lines/) {
                $expected_records = $1;
            }
        }
        $log_file->close;
        $self->throw("The log file $log_path had no information on how many records were expected in the output") unless defined $expected_records;
        
        my $actual_recrods = $out_file->num_records;
        
        if ($actual_recrods == $expected_records) {
            $self->relate_input_to_output($in_path, 'loh_calls', $out_path, { control_sample => $out_file->meta_value('sample_control') });
            
            my $vrtrack              = VRPipe::Schema->create('VRTrack');
            my $graph                = $vrtrack->graph;
            my $output_file_in_graph = $vrtrack->get_file($out_path);
            
            my $md5 = $output_file_in_graph->md5;
            unless ($md5) {
                $md5 = $out_file->file_md5($out_file);
                $output_file_in_graph->md5($md5);
            }
            
            my $fh = $out_file->openr;
            my $sample_source;
            my %calls;
            while (<$fh>) {
                chomp;
                my ($chr, $start, $end, $sample, $count) = split(/\t/, $_);
                my $type = '';
                if ($count =~ /^(\d+)\s+(\S+)/) {
                    $count = $1;
                    $type  = $2;
                    $type =~ s/^\(//;
                    $type =~ s/\)$//;
                }
                
                $sample_source ||= $vrtrack->sample_source($sample);
                my $sample_props = $vrtrack->sample_props_from_string($sample, $sample_source);
                
                push(@{ $calls{ $sample_props->{name} } }, { chr => $chr, start => $start, end => $end, count => $count, type => $type });
            }
            $out_file->close();
            
            my $t = time();
            while (my ($sample, $calls) = each %calls) {
                $calls = [sort { $a->{chr} <=> $b->{chr} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} } @$calls];
                
                $vrtrack->add(
                    'LOH',
                    {
                        md5_sample => $md5 . '_' . $sample,
                        date       => $t,
                        data       => $graph->json_encode($calls)
                    },
                    incoming => [
                        { node_spec => { namespace => 'VRTrack', label => 'Sample', properties => { name => $sample } }, type => 'loh_calls' },
                        { node => $output_file_in_graph, type => 'parsed' }
                    ],
                    enqueue => 1
                );
            }
            $vrtrack->dispatch_queue;
            
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_recrods records were generated in the output file, yet there were $expected_records records expected");
        }
    }
}

1;
