
=head1 NAME

VRPipe::Steps::calculate_bam_spatial_filter - a step

=head1 DESCRIPTION

NB: this step has Sanger-specific hard-coding. It is unlikely to do anything
useful outside the Sanger institute.

Runs spatial_filter to identify regions with spatially correlated errors e.g.
bubbles in aligned BAM files, generating a filter file for subsequent filtering
out or marking. Create a lane-wide filter file from a give tag file, defaults
to tag 168 (the reads aligning to the solexa phiX control) - this allows us to
set consistent thresholds for all bams in the lane. The phix cotrol lane is
assumed to be in irods. Creation of a seperate filter from each bam is possible 
by setting tag_number to zero - probably not to be recommended.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::calculate_bam_spatial_filter with VRPipe::StepRole {
    method options_definition {
        return {
            tag_number => VRPipe::StepOption->create(
                description => 'Generate a lane-wide filter from this tag, eg the tag 168 (phiX control) file. Set to 0 to generate filter from the bam itself',
                optional    => 1, default_value => 168
            ),
            'irods_root' => VRPipe::StepOption->create( description  => 'irods root path under which runParameters xml and Phix Control are held.', optional => 1, default_value => '/seq'),
            samtools_exe               => VRPipe::StepOption->create(description => 'path to your samtools executable',       optional => 1, default_value => 'samtools'),
            samtools_irods_exe         => VRPipe::StepOption->create(description => 'path to your samtools_irods executable', optional => 1, default_value => 'samtools_irods'),
            spatial_filter_exe         => VRPipe::StepOption->create(description => 'path to spatial_filter executable'),
            filter_calibration_options => VRPipe::StepOption->create(
                description   => 'spatial_filter calibration options, excluding tile width and height, and output filter filename',
                optional      => 1,
                default_value => '--region_min_count 122'
            ),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more bam files, grouped by lane metadata if creating lane-wide filter from a specific tag (phiX control)',
                metadata    => {
                    lane        => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    tile_height => 'Tile height',
                    tile_width  => 'Tile width',
                    optional    => ['tile_height', 'tile_width']
                }
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                       = shift;
            my $options                    = $self->options;
            my $samtools                   = $options->{samtools_exe};
            my $samtools_irods_exe         = $options->{samtools_irods_exe};
            my $irods_root                 = $options->{irods_root};
            my $spatial_filter_exe         = $options->{spatial_filter_exe};
            my $tag_number                 = $options->{tag_number};
            my $filter_calibration_options = $options->{filter_calibration_options};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            # Get the lane metadata from the first bam.
            my $bam_file    = $self->inputs->{bam_files}[0];
            my $bam_path    = $bam_file->path;
            my $lane        = $bam_file->metadata->{lane};
            my $tile_width  = $bam_file->metadata->{tile_width} ? $bam_file->metadata->{tile_width} : '0';
            my $tile_height = $bam_file->metadata->{tile_height} ? $bam_file->metadata->{tile_height} : '0';
            
            my $basename = $bam_file->basename;
            $basename .= '.filt';
            
            my $out_file = $self->output_file(output_key => 'filter_files', basename => $basename, type => 'txt', metadata => { lane => $lane });
            my $out_path = $out_file->path;
            
            my $cmd = qq[use VRPipe::Steps::calculate_bam_spatial_filter; VRPipe::Steps::calculate_bam_spatial_filter->generate_filter(source => q[$bam_path], dest => q[$out_path], samtools_exe => '$samtools', samtools_irods_exe => '$samtools_irods_exe', irods_root => '$irods_root', spatial_filter_exe => '$spatial_filter_exe', tag_number => $tag_number, filter_calibration_options => '$filter_calibration_options', lane => '$lane', tile_width => '$tile_width', tile_height => '$tile_height');];
            
            $self->dispatch_vrpipecode($cmd, $req, { output_files => [$out_file] });
        
        };
    }
    
    method outputs_definition {
        return {
            filter_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'filter file identifying bam regions with spatially correlated errors',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs spatial_filter to identify regions with spatially correlated errors e.g. bubbles in aligned BAM files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method generate_filter (ClassName|Object $self: Str|File :$source!, Str|File :$dest!, Str :$samtools_exe!, Str :$samtools_irods_exe!, Str :$irods_root!, Str :$spatial_filter_exe!, Str :$tag_number!, Str :$filter_calibration_options!, Str :$lane!, Num :$tile_width!, Num :$tile_height! ) {
        my $cmd;
        my ($run, undef) = split(/_/, $lane);
        
        if ($tag_number eq '168') {
            # assumption is that the phix files are in irods, eg
            # /seq/8280/8280_4#168.bam
            
            ($lane, undef) = split(/#/, $lane); # strip out tag if present
            $cmd = "$samtools_irods_exe view -u irods:${irods_root}/$run/${lane}#${tag_number}.bam";
        }
        else {
            # read the imported bam directly
            $cmd = "$samtools_exe view -u $source";
        }
        
        unless ($tile_width > 0 && $tile_height > 0) {
            # Tile dimensions are not in the bam metadata, try to get from the RunParameters xml in irods
            my $pipe = "iget ${irods_root}/$run/runParameters.xml - |";
            my $fh;
            open($fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
            while (<$fh>) {
                if (/<TileWidth>(\d+)</) {
                    $tile_width = $1;
                }
                if (/<TileHeight>(\d+)</) {
                    $tile_height = $1;
                }
            }
            close($fh);
            
            $self->throw("Could not iget Tile dimensions from ${irods_root}/$run/RunParameters.xml") unless $tile_width > 0 && $tile_height > 0;
        }
        
        $cmd .= "| $spatial_filter_exe --width $tile_width --height $tile_height $filter_calibration_options -c -F $dest -";
        
        system($cmd) && $self->throw("failed to run [$cmd]");
        
        my $filt_file = VRPipe::File->get(path => $dest);
        $filt_file->update_stats_from_disc;
        unless ($filt_file->s) {
            $self->throw("Failed to generate filter file: '$cmd'");
        }
        return 1;
    }
}

1;
