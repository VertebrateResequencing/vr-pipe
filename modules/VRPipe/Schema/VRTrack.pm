
=head1 NAME

VRPipe::Schema::VRTrack - schemas for tracking biology lab information

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A replacement for our vr-codebase VRTrack mysql-based database schema, this
defines all the things we need to track for our work in Vertebrate Resequencing
at the Sanger. It's pretty-much Sanger-specific.

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

class VRPipe::Schema::VRTrack with VRPipe::SchemaRole {
    method schema_definitions {
        return [
            # general
            {
                label  => 'Group',  # equivalent of old mysql database name, for grouping studies that we will analysis the same way
                unique => [qw(name)]
            },
            {
                label        => 'Study',
                unique       => [qw(id)],
                indexed      => [qw(name accession)],
                keep_history => 1
            },
            {
                label   => 'Taxon',
                unique  => [qw(id)],
                indexed => [qw(common_name)]
            },
            {
                label  => 'Donor',
                unique => [qw(id)]
            },
            {
                label        => 'Sample',
                unique       => [qw(name)],
                indexed      => [qw(id public_name supplier_name accession created_date consent control)],
                keep_history => 1
            },
            {
                label          => 'File',
                unique         => [qw(path)],
                indexed        => [qw(type manual_qc target md5)],
                keep_history   => 1,
                allow_anything => 1
            },
            
            # bams
            {
                label        => 'Library',
                unique       => [qw(id)],
                indexed      => [qw(name tag)],
                keep_history => 1
            },
            {
                label        => 'Lane',
                unique       => [qw(unique)],                             # to be unique but still have correct relationships, this will need to be based on file basename
                required     => [qw(lane)],
                indexed      => [qw(lane run total_reads is_paired_read)],
                keep_history => 1
            },
            {
                label  => 'Alignment',
                unique => [qw(reference)]
            },
            {
                label   => 'EBI_Run',
                unique  => [qw(acc)],
                indexed => [qw(md5)]
            },
            {
                label   => 'EBI_Submission',
                unique  => [qw(acc)],
                indexed => [qw(sub_date)]
            },
            
            # infinium
            {
                label   => 'Beadchip',
                unique  => [qw(id)],
                indexed => [qw(design)]
            },
            {
                label    => 'Section',
                unique   => [qw(unique)], # as per Lane
                required => [qw(section)]
            },
            {
                label  => 'Infinium_Sample',
                unique => [qw(id)]
            },
            {
                label  => 'Infinium_Plate',
                unique => [qw(id)]
            },
            {
                label    => 'Well',
                unique   => [qw(unique)], # as per Lane, but prefix with plate id
                required => [qw(well)]
            },
            
            # analysis files
            {
                label  => 'Analysis',
                unique => [qw(uuid)]
            },
            {
                label    => 'Collection',
                unique   => [qw(path)],
                required => [qw(date)]
            },
            
            # QC bam stats
            {
                label    => 'Bam_Stats',
                unique   => [qw(uuid)],
                required => [qw(reads bases)],
                optional => [qw(mean_coverage)]
            },
        ];
    }
}

1;
