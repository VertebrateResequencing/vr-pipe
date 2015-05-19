
=head1 NAME

VRPipe::Steps::vrtrack_populate_from_graph_db - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::vrtrack_populate_from_graph_db extends VRPipe::Steps::vrtrack_populate_from_vrpipe_metadata {
    use VRPipe::Schema;
    use VRPipe::Steps::vrtrack_update_mapstats;
    
    method body_sub {
        return sub {
            my $self        = shift;
            my $opts        = $self->options;
            my $db          = $opts->{vrtrack_db};
            my $storage_dir = $opts->{vrlane_storage_dir};
            my $req         = $self->new_requirements(memory => 500, time => 1);
            foreach my $file (@{ $self->inputs->{files} }) {
                my $fid = $file->id;
                my $cmd = "use VRPipe::Steps::vrtrack_populate_from_graph_db; VRPipe::Steps::vrtrack_populate_from_graph_db->populate_lane_from_file(db => q[$db], file => $fid";
                if ($storage_dir) {
                    $cmd .= ", storage_dir => q[$storage_dir]";
                }
                $cmd .= ");";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method description {
        return "Populate or update VRTrack lanes and mapstats based on information the irods all_with_warehouse_metadata datasource, npg_cram_stats_parser step and plot_bamstats step put in the graph database.";
    }
    
    method populate_lane_from_file (ClassName|Object $self: Str :$db!, Int :$file!, Str|Dir :$storage_dir?) {
        $file = VRPipe::File->get(id => $file);
        
        my $schema = VRPipe::Schema->create("VRPipe");
        my $graph_file = $schema->get('File', { path => $file->protocolless_path, protocol => $file->protocol });
        $self->throw($file->path . " was not in the graph db") unless $graph_file;
        
        # create/update lane hierarchy
        my $props = $graph_file->properties(flatten_parents => 1);
        my $lane = $props->{vrtrack_lane_unique} || $self->throw($graph_file->path . " had no associated lane");
        
        my %conversion = (
            total_reads             => 'lane_total_reads',
            study_id                => 'study_id',
            is_paired_read          => 'lane_is_paired_read',
            manual_qc               => 'manual_qc',
            md5                     => 'md5',
            library                 => 'library_name',
            library_id              => 'library_id',
            infinium_sample         => 'infinium_sample_id',
            sample                  => 'sample_name',
            beadchip                => 'beadchip_id',
            beadchip_section        => 'section_section',     #*** or section_unique?
            sample_id               => 'sample_id',
            study_title             => 'study_name',
            study_accession_number  => 'study_accession',
            sample_supplier_name    => 'sample_supplier_name',
            sample_control          => 'sample_control',
            taxon_id                => 'taxon_id',
            sample_common_name      => 'taxon_common_name',
            sample_cohort           => 'sample_cohort',
            sample_accession_number => 'sample_accession',
            public_name             => 'sample_public_name'
        );
        
        my $meta = {};
        while (my ($key, $graph_key) = each %conversion) {
            my $val = $props->{ 'vrtrack_' . $graph_key };
            next unless defined $val;
            $meta->{$key} = $val;
        }
        
        if ($storage_dir) {
            #*** analysis stuff should also be queriable from the graph db, but
            # I just cheat for now, taking it from the file metadata that the
            # irods datasource adds
            my $file_meta = $file->metadata;
            if (defined $file_meta->{analysis_uuid}) {
                $meta->{analysis_uuid}        = $file_meta->{analysis_uuid};
                $meta->{irods_analysis_files} = $file_meta->{irods_analysis_files};
            }
        }
        $self->meta_to_lane(db => $db, file => $file, meta => $meta, lane => $lane, $storage_dir ? (storage_dir => $storage_dir) : ());
        
        # create/update mapstats and graphs
        my %plot_files;
        my ($bs) = sort { $b->date <=> $a->date } $graph_file->related(outgoing => { namespace => "VRTrack", label => "Bam_Stats", max_depth => 4 });
        $self->throw("No Bam_Stats node found related to " . $graph_file->path) unless $bs;
        my ($bs_file) = $bs->related(incoming => { type => "summary_stats" });
        $self->throw("No stats file node found related to Bam_Stats node for " . $graph_file->path) unless $bs_file;
        foreach my $pf ($bs_file->related(outgoing => { type => "bamstats_plot" })) {
            $plot_files{ $pf->path } = $pf->property('caption');
        }
        
        my $mode            = $bs->mode;
        my $meta_key_prefix = $mode ne 'normal' ? $mode . '_' : '';
        my $targets_mode    = $mode eq 'targeted' ? 1 : 0;
        
        %conversion = (
            filtered_reads     => 'filtered sequences',
            reads              => 'sequences',                     #*** 'raw total sequences' is total records, 'sequences' is the 0x900 count, but will this be true in samtools 1.3+?
            bases              => 'total length',
            reads_mapped       => 'reads mapped',
            reads_paired       => 'reads paired',
            bases_mapped_c     => 'bases mapped (cigar)',
            error_rate         => 'error rate',
            rmdup_reads_mapped => 'reads mapped after rmdup',
            rmdup_bases_mapped => 'bases mapped after rmdup',
            bases_trimmed      => 'bases trimmed',
            mean_insert_size   => 'insert size average',
            sd_insert_size     => 'insert size standard deviation',
            avg_read_length    => 'average length'
        );
        if ($targets_mode) {
            foreach my $cov (1, 2, 5, 10, 20, 50, 100) {
                $conversion{"targeted_bases_of_${cov}X_coverage"} = "bases of ${cov}X coverage";
            }
        }
        
        $meta = {};
        while (my ($key, $graph_key) = each %conversion) {
            my $prefixed_key = $meta_key_prefix . $key;
            my $val          = $bs->property($graph_key);
            next unless defined $val;
            $meta->{$prefixed_key} = $val;
        }
        $meta->{paired}        = $bs->property('reads properly paired') ? 1 : 0;
        $meta->{npg_qc_status} = $graph_file->property('manual_qc')     ? 1 : 0;
        
        VRPipe::Steps::vrtrack_update_mapstats->meta_to_mapstats(meta => $meta, plot_files => \%plot_files, lane_file => $graph_file, db => $db, lane => $lane, meta_key_prefix => $meta_key_prefix, targets_mode => $targets_mode);
    }
}

1;
