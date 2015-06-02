
=head1 NAME

VRPipe::Steps::vcf_genotype_comparison - a step

=head1 DESCRIPTION

Compare the genotypes of the samples in a VCF using bcftools gtcheck. Adds
'genotype_maximum_deviation' metadata to the VCF file indicating if all its
samples were similar to each other. And stores the details of all pairwise
discordance results in the graph database.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2015 Genome Research Limited.

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

class VRPipe::Steps::vcf_genotype_comparison extends VRPipe::Steps::bcftools {
    use VRPipe::Schema;
    
    around options_definition {
        return {
            %{ $self->$orig },
            bcftools_gtcheck_options => VRPipe::StepOption->create(
                description => 'options to bcftools gtcheck (-g option is not valid for this step)',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'vcf files; genotype comparison will be done on each independently'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self         = shift;
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            my $gtcheck_opts = $options->{bcftools_gtcheck_options};
            if ($gtcheck_opts =~ /gtcheck| -g| --genotypes/) {
                $self->throw("bcftools_gtcheck_options should not include the gtcheck subcommand or the -g option");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => $self->bcftools_version_string,
                    summary => "bcftools gtcheck $gtcheck_opts \$vcf_file > \$gtcheck_file.gtypex"
                )
            );
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my $vcf_path = $vcf->path;
                my $meta     = $vcf->metadata;
                
                my $gtypex_file = $self->output_file(
                    output_key => 'bcftools_gtcheck_files',
                    basename   => $vcf->basename . '.gtypex',
                    type       => 'txt',
                    metadata   => $vcf->metadata
                );
                my $gtypex_path = $gtypex_file->path;
                my $cmd         = qq[$bcftools_exe gtcheck $gtcheck_opts $vcf_path > $gtypex_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_genotype_comparison', 'compare_genotypes', [$cmd, $req, { output_files => [$gtypex_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bcftools_gtcheck_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'file of genotype concurrence scores calculated by bcftools gtcheck',
                metadata    => { genotype_maximum_deviation => "maximum deviation in genotype and the sample causing that deviation" }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Compare the genotypes of the samples in a VCF to each other to confirm they come from the same individual using bcftools gtcheck";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method compare_genotypes (ClassName|Object $self: Str $cmd_line) {
        my ($vcf_path, $output_path) = $cmd_line =~ /(\S+) > (\S+)$/;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $vcf_file    = VRPipe::File->get(path => $vcf_path);
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $vrtrack              = VRPipe::Schema->create('VRTrack');
        my $output_file_in_graph = $vrtrack->add_file($output_path);
        $self->relate_input_to_output($vcf_path, 'genotypes_compared', $output_file_in_graph);
        
        my $md5 = $output_file_in_graph->md5;
        unless ($md5) {
            $md5 = $output_file->file_md5($output_file);
            $output_file_in_graph->md5($md5);
        }
        
        my $graph = $vrtrack->graph;
        
        my $fh    = $output_file->openr;
        my $count = 0;
        my %pairs;
        my $sample_source;
        my $t      = time();
        my $max_cn = 0;
        my $gmd;
        while (<$fh>) {
            if (/^SM\s/) {
                # gmd used to be parsable direct from MD line, but now not all
                # versions of bcftools produce it; they all produce SMs, so we
                # just calculate it from here
                my (undef, $discordance, undef, undef, $sample) = split;
                
                if ($discordance >= $max_cn) {
                    $max_cn = $discordance;
                    $gmd    = "$discordance:$sample";
                }
            }
            elsif (/^CN\s/) {
                chomp;
                my (undef, $discordance, $num_of_sites, $avg_min_depth, $sample_i, $sample_j) = split;
                
                # there could be multiple rows with the same pair of samples and
                # we need a unique identifier for each one; in the file they are
                # uniqufied by adding some number prefix to one of the samples,
                # but we don't like that
                $sample_i =~ s/^\d+\://;
                $sample_j =~ s/^\d+\://;
                my $prop_key = join('.', sort ($sample_i, $sample_j));
                my $i = 0;
                while (1) {
                    $i++;
                    my $this_key = $prop_key . '.' . $i;
                    next if exists $pairs{$this_key};
                    $pairs{$this_key} = 1;
                    $prop_key = $this_key;
                    last;
                }
                
                # before we can store the result we need to find the sample
                # nodes in the graph database under the VRTrack schema;
                # complication is that the sample names in the file could be
                # name, public_name or some combination of both (or indeed
                # anything else). Further problem is that since this is a
                # wrapped cmd, we can't even pass in the info of what was used
                # as an argument... for now we try out the obvious possibilities
                # until found
                unless (defined $sample_source) {
                    my $sample_node;
                    if ($sample_i =~ /^(.+)_([^_]+)$/) {
                        # public_name+sample
                        $sample_node = $vrtrack->get('Sample', { public_name => $1, name => $2 });
                        if ($sample_node) {
                            $sample_source = 'public_name+sample';
                        }
                    }
                    if (!$sample_node) {
                        # sample
                        $sample_node = $vrtrack->get('Sample', { name => $sample_i });
                        if ($sample_node) {
                            $sample_source = 'sample';
                        }
                    }
                    if (!$sample_node) {
                        # ... give up for now
                        $self->throw("Couldn't find a Sample node for $sample_i in the graph database");
                    }
                }
                
                my @props;
                foreach my $identifer ($sample_i, $sample_j) {
                    my $sample_props;
                    if ($sample_source eq 'public_name+sample') {
                        my ($public_name, $sample) = $identifer =~ /^(.+)_([^_]+)$/;
                        $sample_props = { public_name => $public_name, name => $sample };
                    }
                    elsif ($sample_source eq 'sample') {
                        $sample_props = { name => $identifer };
                    }
                    push(@props, $sample_props);
                }
                
                foreach my $i (1 .. 2) {
                    # we add in enqueue mode which doesn't add straight away,
                    # but only when we call dispatch_queue() every 10000 loops -
                    # for database efficiency. The first time will create a
                    # discordance node for this sample; subsequent times will
                    # just add a new property keyed on $prop_key with the
                    # discordance results in an array, with the other sample
                    # name in there as well to make it easy to know which
                    # results are vs which sample without additional db queries
                    # later
                    $vrtrack->add(
                        'Discordance',
                        {
                            # we unique on md5 of the gtypex file and sample
                            # name, so don't change the graph (other than date)
                            # if we've parsed this file before, and don't
                            # alter results on a node for a different gtypex
                            # file or sample
                            md5_sample => $md5 . '_' . $props[$i - 1]->{name},
                            date       => $t,
                            type       => $num_of_sites < 100 ? 'fluidigm' : 'genotype',
                            $prop_key => ["$discordance", "$num_of_sites", "$avg_min_depth", $i == 1 ? $props[1]->{name} : $props[0]->{name}]
                        },
                        incoming => [
                            { node_spec => { namespace => 'VRTrack', label => 'Sample', properties => $props[$i - 1] }, type => 'discordance' },
                            { node => $output_file_in_graph, type => 'parsed' }
                        ],
                        enqueue => 1
                    );
                }
                
                $count++;
                if ($count == 10000) {
                    $vrtrack->dispatch_queue;
                    $count = 0;
                }
            }
        }
        $output_file->close;
        $vrtrack->dispatch_queue;
        
        if ($gmd) {
            foreach my $file ($vcf_file, $output_file) {
                $file->add_metadata({ genotype_maximum_deviation => $gmd });
            }
        }
    }
}

1;
