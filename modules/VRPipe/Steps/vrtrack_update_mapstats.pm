
=head1 NAME

VRPipe::Steps::vrtrack_update_mapstats - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::vrtrack_update_mapstats extends VRPipe::Steps::vrtrack_update {
    around options_definition {
        return {
            %{ $self->$orig },
            exome_targets_file => VRPipe::StepOption->create(
                description => 'absolute path to a file describing the targets/baits used for exome pulldown (tab-delimited [chr,start,end], where start is 1-based, and end is inclusive)',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                description => 'bam file with associated bamcheck statistics in the metadata',
                max_files   => -1,
                metadata    => { lane => 'lane name (a unique identifer for this sequencing run, aka read group)' }
            ),
            bamcheck_plots => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'png files produced by plot-bamcheck, with a caption in the metadata',
                min_files   => 11,
                max_files   => -1,
                metadata    => {
                    source_bam => 'the bam file this plot was made from',
                    caption    => 'the caption of this plot'
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $opts          = $self->options;
            my $db            = $opts->{vrtrack_db};
            my $targeted_mode = $opts->{exome_targets_file} ? 1 : 0;
            my $req           = $self->new_requirements(memory => 500, time => 1);
            
            my %bam_plots;
            foreach my $plot_file (@{ $self->inputs->{bamcheck_plots} }) {
                my $source_bam = $plot_file->metadata->{source_bam};
                $bam_plots{$source_bam}->{dir} = $plot_file->dir;
                push(@{ $bam_plots{$source_bam}->{files} }, $plot_file->basename);
            }
            
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $bam_path        = $bam_file->path;
                my $lane            = $bam_file->metadata->{lane};
                my $these_bam_plots = $bam_plots{$bam_path};
                my $cmd             = "use VRPipe::Steps::vrtrack_update_mapstats; VRPipe::Steps::vrtrack_update_mapstats->update_mapstats(db => q[$db], bam => q[$bam_path], lane => q[$lane], targets_mode => $targeted_mode, plot_dir => q[$these_bam_plots->{dir}], plots => [qw[@{$these_bam_plots->{files}}]]);";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method description {
        return "Add the bamcheck QC statistics and graphs to the VRTrack database, so that they're accessible with QCGrind etc.";
    }
    
    method update_mapstats (ClassName|Object $self: Str :$db!, Str|File :$bam!, Str :$lane!, Str|Dir :$plot_dir!, ArrayRef :$plots!, Bool :$targets_mode?) {
        my $bam_file = VRPipe::File->get(path => $bam);
        my $meta = $bam_file->metadata;
        my %plot_files;
        foreach my $plot_basename (@$plots) {
            my $file = VRPipe::File->get(path => file($plot_dir, $plot_basename));
            $plot_files{ $file->path } = $file->metadata->{caption};
        }
        $bam_file->disconnect;
        
        # do we need to sanity check the .dict file vs the SQ lines in the bam
        # header?
        
        # get the lane and mapstats object from VRTrack
        my $vrtrack = $self->get_vrtrack(db => $db);
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("No lane named '$lane' in database '$db'");
        my $mapstats = $vrlane->latest_mapping;
        
        my $worked = $vrtrack->transaction(
            sub {
                unless ($mapstats) {
                    $mapstats = $vrlane->add_mapping();
                    #*** do we need to fill in mapper and assembly?
                }
                
                # fill in the mapstats based on our bam metadata
                my $meta_key_prefix = $targets_mode ? 'targeted_' : '';
                
                $mapstats->raw_reads($meta->{ $meta_key_prefix . 'filtered_reads' });
                my $raw_bases = $meta->{ $meta_key_prefix . 'bases' };
                $mapstats->raw_bases($raw_bases);
                $mapstats->reads_mapped($meta->{ $meta_key_prefix . 'reads_mapped' });
                $mapstats->reads_paired($meta->{ $meta_key_prefix . 'reads_paired' });
                $mapstats->bases_mapped($meta->{ $meta_key_prefix . 'bases_mapped_c' });
                $mapstats->error_rate($meta->{ $meta_key_prefix . 'error_rate' });
                $mapstats->rmdup_reads_mapped($meta->{ $meta_key_prefix . 'rmdup_reads_mapped' });
                $mapstats->rmdup_bases_mapped($meta->{ $meta_key_prefix . 'rmdup_bases_mapped' }); # ideally this would be rmdup_bases_mapped_c, but we no longer calculate this with bamcheck -d
                $mapstats->clip_bases($raw_bases - $meta->{ $meta_key_prefix . 'bases_trimmed' });
                $mapstats->mean_insert($meta->{ $meta_key_prefix . 'mean_insert_size' });
                $mapstats->sd_insert($meta->{ $meta_key_prefix . 'sd_insert_size' });
                
                if ($targets_mode) {
                    foreach my $cov (1, 2, 5, 10, 20, 50, 100) {
                        my $these_bases = $meta->{"targeted_bases_of_${cov}X_coverage"} || next;
                        my $method = "target_bases_${cov}X";
                        $mapstats->$method(sprintf("%0.2f", $these_bases / $raw_bases));
                    }
                    $mapstats->mean_target_coverage($meta->{targeted_mean_coverage}) if $meta->{targeted_mean_coverage};
                    $mapstats->target_bases_mapped($meta->{targeted_bases_mapped_c});
                }
                
                # add the images
                while (my ($path, $caption) = each %plot_files) {
                    my $img = $mapstats->add_image_by_filename($path);
                    $img->caption($caption);
                    $img->update;
                }
                
                $mapstats->update;
                
                # say that the file is imported
                my $vrfile = $vrlane->get_file_by_name($bam_file->basename);
                if ($vrfile) {
                    $vrfile->is_processed(import => 1);
                    $vrfile->md5($bam_file->md5);
                    $vrfile->update;
                    $vrfile->is_processed(mapped => 1);
                    $vrfile->update;
                }
                
                # also update the lane, only changing qc_status if it is at the
                # default status of no_qc
                $vrlane->raw_bases($meta->{bases});
                $vrlane->raw_reads($meta->{reads});
                $vrlane->is_paired($meta->{paired} ? 1 : 0);
                $vrlane->read_len(int($meta->{avg_read_length}));
                if ($vrlane->qc_status eq 'no_qc') {
                    $vrlane->qc_status('pending');
                }
                $vrlane->is_processed(import => 1);
                $vrlane->update;
                $vrlane->is_processed(mapped => 1);
                $vrlane->update;
                $vrlane->is_processed(qc => 1);
                $vrlane->update;
            }
        );
        
        unless ($worked) {
            $self->throw($vrtrack->{transaction_error});
        }
    }
}

1;
