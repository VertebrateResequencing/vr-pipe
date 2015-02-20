
=head1 NAME

VRPipe::Steps::samtools_bam_stats - a step

=head1 DESCRIPTION

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

class VRPipe::Steps::samtools_bam_stats with VRPipe::StepRole {
    use VRPipe::Parser;
    use VRPipe::Schema;
    
    method options_definition {
        return {
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            ),
            samtools_stats_options => VRPipe::StepOption->create(
                description => 'options to samtools stats, excluding -r and -t (which are set by reference_fasta and exome_targets_file options)',
                optional    => 1
            ),
            reference_fasta    => VRPipe::StepOption->create(description => 'absolute path to genome reference file'),
            exome_targets_file => VRPipe::StepOption->create(
                description => 'absolute path to a file describing the targets/baits used for exome pulldown (tab-delimited [chr,start,end], where start is 1-based, and end is inclusive)',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', description => 'bam files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            my $samtools_exe = $options->{samtools_exe};
            my $opts         = VRPipe::Steps::samtools_bam_stats->get_stats_options($options);
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $ifile      = $bam_file->path;
                my $stats_file = $self->output_file(output_key => 'stats_files', basename => $ifile->basename . '.bamstats', type => 'txt');
                my $ofile      = $stats_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::samtools_bam_stats', 'stats', ["$samtools_exe stats $opts $ifile > $ofile", $req, { output_files => [$stats_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            stats_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'the output of samtools stats on a bam',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Creates a file of stats on a bam file, and also stores the summary stats in the graph database under the VRTrack schema.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method get_stats_options (ClassName|Object $self: HashRef $options!) {
        my $ref = file($options->{reference_fasta});
        $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
        my $opts = '-r ' . $ref;
        
        my $targets = $options->{exome_targets_file};
        if ($targets) {
            $targets = file($targets);
            $self->throw("exome_targets_file must be an absolute path") unless $targets->is_absolute;
            $opts .= ' -t ' . $targets;
        }
        
        my $user_opts = $options->{samtools_stats_options};
        if ($user_opts) {
            if ($user_opts =~ /-r|-t/) {
                $self->throw("neither -r nor -t should be supplied as samtools stats options");
            }
            
            $opts .= ' ' . $user_opts;
        }
        
        return $opts;
    }
    
    method stats (ClassName|Object $self: Str $cmd_line) {
        my ($opts, $bam_path, $stats_path) = $cmd_line =~ /^\S+ stats (.*) (\S+) > (\S+)$/;
        $opts ||= '';
        $bam_path || $self->throw("bad cmd line [$cmd_line]");
        my $bam_file   = VRPipe::File->get(path => $bam_path);
        my $stats_file = VRPipe::File->get(path => $stats_path);
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $stats_file->update_stats_from_disc(retries => 3);
        if ($stats_file->s) {
            # parse the stats file with the old bamcheck parser
            my $parser = VRPipe::Parser->create('bamcheck', { file => $stats_file });
            $stats_file->disconnect;
            
            # we'll store results in the graph database under the VRTrack schema
            my $vrtrack = VRPipe::Schema->create('VRTrack');
            
            # before adding graph node for the stats file, make sure we have
            # nodes for the whole hierarchy
            $parser = VRPipe::Parser->create('bam', { file => $bam_file });
            my %rg_info = $parser->readgroup_info();
            my @rgs     = keys %rg_info;
            my @lanes;
            foreach my $i (0 .. $#rgs) {
                my $info              = $rg_info{ $rgs[$i] };
                my $unique_lane_id    = $info->{PU} || $rgs[$i];
                my $changed_lane_name = 0;
                if ("$unique_lane_id" eq "1") {
                    # call the lane something we can be most sure is maximally
                    # unique
                    my $md5 = $bam_file->md5;
                    unless ($md5) {
                        $bam_file->update_md5;
                        $md5 = $bam_file->md5;
                    }
                    $unique_lane_id    = $md5 . $i;
                    $changed_lane_name = 1;
                }
                
                my $nodes = $vrtrack->ensure_sequencing_hierarchy(
                    lane    => $unique_lane_id,
                    library => $info->{LB},
                    sample  => $info->{SM},
                    study   => $info->{DS}
                );
                if ($nodes->{library}) {
                    $nodes->{library}->add_properties({ center_name => $info->{CN}, platform => $info->{PL} });
                }
                if ($nodes->{study}) {
                    $vrtrack->add('Group', { name => 'all_studies' }, outgoing => { type => 'has', node => $nodes->{study} });
                }
                $nodes->{lane}->lane($info->{PU} || $rgs[$i]) if $changed_lane_name;
                push(@lanes, $nodes->{lane});
            }
            
            # now add the stats file node and relate it with what we can
            my $vrsource = $vrtrack->get_file($bam_path);
            my $vrfile   = $vrtrack->add_file($stats_path);
            if ($vrsource) {
                $vrsource->relate_to($vrfile, 'parsed');
            }
            elsif (@lanes) {
                foreach my $vrlane (@lanes) {
                    $vrlane->relate_to($vrfile, 'stats');
                }
            }
            
            if ($vrfile) {
                $self->result_nodes([$vrfile]);
                
                # parse out stats from the stats file (don't use our bamcheck
                # parser so we just store whatever new or changed field is in
                # the stats file)
                my $fh    = $stats_file->openr;
                my %stats = ();
                while (<$fh>) {
                    if (/^SN\s+([^:]+):\s+(\S+)/) {
                        $stats{$1} = $2;
                        
                        # store reads metadata on the file for compatibility
                        # with old steps
                        if ($1 eq 'raw total sequences') {
                            $bam_file->add_metadata({ reads => $2 });
                        }
                    }
                }
                $stats_file->close;
                
                my $mode = 'normal';
                if ($cmd_line =~ /-d/) {
                    $mode = 'rmdup';
                }
                elsif ($cmd_line =~ /-t/) {
                    $mode = 'targeted';
                }
                
                unless ($mode eq 'rmdup') {
                    # when not using -d we can still add some rmdup stats
                    my $reads_duplicated = $stats{'reads duplicated'};
                    $stats{'reads after rmdup'}        = $stats{'raw total sequences'} - $reads_duplicated;
                    $stats{'reads mapped after rmdup'} = $stats{'reads mapped'} - $reads_duplicated;
                    my $bases_duplicated = $stats{'bases duplicated'};
                    $stats{'bases after rmdup'} = $stats{'total length'} - $bases_duplicated;
                    
                    # Estimate rmdup_bases_mapped CIGAR (excluding soft clipped
                    # and for exomes parts of the reads which fall outside of
                    # the target regions ) Calculate fraction
                    # bases_mapped/bases_duplicated and apply this to
                    # bases_mapped_cigar
                    my $bases_mapped_cigar = $stats{'bases mapped (cigar)'};
                    my $bases_dup_fraction = $stats{'bases mapped'} ? $bases_duplicated / $stats{'bases mapped'} : 0;
                    my $bases_dup_estimate = int($bases_mapped_cigar * $bases_dup_fraction);
                    $stats{'bases mapped after rmdup'} = $bases_mapped_cigar - $bases_dup_estimate;
                    
                    # coverage stats we just don't bother having if in -d mode;
                    # we use the bamcheck parser for this since it's non-trivial
                    $parser = VRPipe::Parser->create('bamcheck', { file => $stats_file });
                    $stats{'mean coverage'} = $parser->mean_coverage;
                    foreach my $cov (1, 2, 5, 10, 20, 50, 100) {
                        $stats{"bases of ${cov}X coverage"} = $parser->cumulative_coverage($cov);
                    }
                }
                
                # store the parsed stats in a Bam_Stats node attached to the
                # file; because they are given an automatic new uuid as the
                # primary identifier for the node, we always end up adding a
                # new bam_stats node instead of updating an existing one. So we
                # store the current date so the user can pick the most recent
                # one
                my $bs = $vrtrack->add(
                    'Bam_Stats',
                    {
                        mode    => $mode,
                        options => $opts,
                        date    => time(),
                        %stats
                    },
                    incoming => { type => 'summary_stats', node => $vrfile }
                );
            }
        }
        else {
            $self->throw("$stats_path failed to be made");
        }
        return 1;
    }
}

1;
