
=head1 NAME

VRPipe::Steps::combine_bcftools_cnvs - a step

=head1 DESCRIPTION

This step runs cmp-cnvs.pl to combine the summary results from bcftools cnv 
across the samples.

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

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

class VRPipe::Steps::combine_bcftools_cnvs with VRPipe::StepRole {
    use VRPipe::Schema;
    
    method options_definition {
        return {
            compare_cnv_script => VRPipe::StepOption->create(
                description   => 'path to your cmp-cnvs.pl which combines bcftools outputs - skips if not provided',
                optional      => 1,
                default_value => 'cmp-cnvs.pl',
            ),
            compare_cnv_options => VRPipe::StepOption->create(
                description   => 'options to cmp-cnvs.pl which combines output from bcftools',
                optional      => 1,
                default_value => '-q 2'
            ),
        };
    }
    
    method inputs_definition {
        return {
            summary_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'summary files from bcftools cnv',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $cmp_cnvs_exe  = $options->{compare_cnv_script};
            my $cmp_cnvs_opts = $options->{compare_cnv_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'cmp-cnvs.pl',
                    version => 0,
                    summary => "cmp-cnvs.pl $cmp_cnvs_opts \@dirs > combined_cnvs.txt"
                )
            );
            
            my $merged_meta   = $self->combined_metadata($self->inputs->{summary_files});
            my $cmp_file      = $self->output_file(output_key => 'cmp_file', basename => 'combined_cnvs.txt', type => 'txt', metadata => $merged_meta);
            my $cmp_file_path = $cmp_file->path->stringify;
            
            my @dirs;
            my @input_strings;
            foreach my $file (@{ $self->inputs->{summary_files} }) {
                push(@dirs,          $file->dir);
                push(@input_strings, $file->path->stringify);
            }
            $self->relate_input_to_output(\@input_strings, 'combined', $cmp_file_path);
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            my $cmd_line = "$cmp_cnvs_exe $cmp_cnvs_opts @dirs > $cmp_file_path";
            $self->dispatch_wrapped_cmd('VRPipe::Steps::combine_bcftools_cnvs', 'combine_cnvs', [$cmd_line, $req, { output_files => [$cmp_file] }]);
        };
    }
    
    method outputs_definition {
        return {
            cmp_file => VRPipe::StepIODefinition->create(
                type            => 'txt',
                min_files       => 1,
                max_files       => 1,
                description     => 'combined results from cmp-cnvs.pl',
                check_existence => 0,
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Combine bcftools cnv summary results across the samples.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method combine_cnvs (ClassName|Object $self: Str $cmd_line) {
        my ($output_path) = $cmd_line =~ /\S > (\S+)$/;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $vrtrack              = VRPipe::Schema->create('VRTrack');
        my $output_file_in_graph = $vrtrack->get_file($output_path);
        
        my $md5 = $output_file_in_graph->md5;
        unless ($md5) {
            $md5 = $output_file->file_md5($output_file);
            $output_file_in_graph->md5($md5);
        }
        
        my $graph = $vrtrack->graph;
        
        my $fh = $output_file->openr;
        my (@samples, %results, %node_specs);
        while (<$fh>) {
            next if /^#/;
            chomp;
            my @cols = split;
            if ($cols[0] eq 'SM') {
                my $sample_source = $vrtrack->sample_source($cols[1]);
                foreach my $sample_str (@cols[1 .. $#cols]) {
                    my $sample_props = $vrtrack->sample_props_from_string($sample_str, $sample_source);
                    push(@samples, $sample_props->{name});
                    $node_specs{ $sample_props->{name} } = { namespace => 'VRTrack', label => 'Sample', properties => $sample_props };
                }
            }
            elsif ($cols[0] eq 'RG') {
                my $chr = $cols[1];
                foreach my $i (7 .. $#cols) {
                    push(@{ $results{ $samples[$i - 6] }->{ $cols[0] } }, { chr => $cols[1], start => $cols[2], end => $cols[3], length => $cols[4], quality => $cols[5], cn => $cols[$i] });
                }
            }
            else {
                foreach my $i (1 .. $#cols) {
                    $results{ $samples[$i - 1] }->{ $cols[0] } = $cols[$i];
                }
            }
        }
        $output_file->close();
        
        my $t = time();
        while (my ($sample, $results) = each %results) {
            $vrtrack->add(
                'CNVs',
                {
                    md5_sample => $md5 . '_' . $sample,
                    date       => $t,
                    data       => $graph->json_encode($results)
                },
                incoming => [
                    { node_spec => $node_specs{$sample},  type => 'cnv_calls' },
                    { node      => $output_file_in_graph, type => 'parsed' }
                ],
                enqueue => 1
            );
        }
        $vrtrack->dispatch_queue;
    }
}

1;
