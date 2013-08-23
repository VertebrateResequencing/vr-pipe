
=head1 NAME

VRPipe::Steps::pluritest_vrtrack_update_images - a step

=head1 DESCRIPTION

Inserts the pluritest images into the images table of the vrtrack database. The
'tag' id for these plots is stored in the mapstats table as reads_mapped.

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::pluritest_vrtrack_update_images extends VRPipe::Steps::vrtrack_update {
    around options_definition {
        return {
            %{ $self->$orig },
        };
    }
    
    method inputs_definition {
        return {
            pluritest_plots => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'png files produced by pluritest R script to assist in determination of pluripotency of stem cell lines',
                min_files   => 5,
                max_files   => -1,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            my $opts = $self->options;
            my $db   = $opts->{vrtrack_db};
            my $req  = $self->new_requirements(memory => 500, time => 1);
            
            my %plots;
            my $analysis_tag;
            foreach my $plot_file (@{ $self->inputs->{pluritest_plots} }) {
                my $lanes = $plot_file->metadata->{lanes};
                $analysis_tag = $plot_file->metadata->{merge_tag_id};
                $plots{$lanes}->{dir} = $plot_file->dir;
                push(@{ $plots{$lanes}->{files} }, $plot_file->basename);
            }
            
            foreach my $lanes (keys %plots) {
                my @analysis_lanes = split(',', $lanes);
                for my $lane (@analysis_lanes) {
                    my $these_bam_plots = $plots{$lanes};
                    my $cmd             = "use VRPipe::Steps::pluritest_vrtrack_update_images; VRPipe::Steps::pluritest_vrtrack_update_images->update_images(db => q[$db], lane => q[$lane], tag => q[$analysis_tag], plot_dir => q[$these_bam_plots->{dir}], plots => [qw[@{$these_bam_plots->{files}}]]);";
                    $self->dispatch_vrpipecode($cmd, $req);
                }
            }
        };
    }
    
    method outputs_definition {
        return {};
    }
    
    method description {
        return "Inserts the pluritest images into the images table of the vrtrack database. The 'tag' id for these plots is stored in the mapstats table as reads_mapped.";
    }
    
    method update_images (ClassName|Object $self: Str :$db!, Str :$lane!, Str :$tag!, Str|Dir :$plot_dir!, ArrayRef :$plots!) {
        my %plot_files;
        foreach my $plot_basename (@$plots) {
            my $file = VRPipe::File->get(path => file($plot_dir, $plot_basename));
            $plot_files{ $file->path } = $tag;
        }
        
        # get the lane and mapstats object from VRTrack
        my $vrtrack = $self->get_vrtrack(db => $db);
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("No lane named '$lane' in database '$db'");
        my $mapstats = $vrlane->latest_mapping;
        
        my $worked = $vrtrack->transaction(
            sub {
                unless ($mapstats) {
                    $mapstats = $vrlane->add_mapping();
                }
                
                # add the images
                while (my ($path, $caption) = each %plot_files) {
                    my $img = $mapstats->add_image_by_filename($path);
                    $img->caption($caption);
                    $img->update;
                    $mapstats->genotype_expected($caption);
                }
                $mapstats->update;
                
                $vrlane->is_processed(import => 1);
                $vrlane->update;
                $vrlane->is_processed(mapped => 1);
                $vrlane->update;
            }
        );
        
        unless ($worked) {
            $self->throw($vrtrack->{transaction_error});
        }
    }
}

1;
