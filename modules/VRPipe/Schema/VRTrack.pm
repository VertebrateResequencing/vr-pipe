
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
    my $vrpipe_schema;
    
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
                label    => 'Donor',
                unique   => [qw(id)],
                optional => [qw(example_sample last_sample_added_date)] # psuedo fields not stored in db, but you could calculate these when retrieving a donor
            },
            {
                label        => 'Sample',
                unique       => [qw(name)],
                indexed      => [qw(id public_name supplier_name accession created_date consent control)],
                keep_history => 1
            },
            
            # bams
            {
                label        => 'Library',
                unique       => [qw(id)],
                indexed      => [qw(name tag platform center_name)],
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
                label          => 'Bam_Stats',
                unique         => [qw(uuid)],
                required       => [qw(mode options), 'raw total sequences'],
                allow_anything => 1
            },
            
            # gender (for expected or calculated from sequenom/fluidgm)
            {
                label   => 'Gender',
                unique  => [qw(source_gender_md5)],
                indexed => [qw(source gender)]
            },
            
            # genotype comparison results
            {
                label    => 'Discordance',
                unique   => [qw(sample_pair)],
                required => [qw(discordance num_of_sites avg_min_depth)]
            },
        ];
    }
    
    # without enforce => 1, if the hierarchy already exists but is different,
    # no changes are made to the database and only nodes up to the difference
    # are returned
    method ensure_sequencing_hierarchy (Str :$lane!, Str :$library!, Str :$sample!, Bool :$enforce = 0, Str :$study?, Str :$taxon?, Str :$group?) {
        # get and only then try adding Lane, because we don't want to override
        # an existing lane property
        my $vrlane = $self->get('Lane', { unique => $lane });
        $vrlane ||= $self->add('Lane', { unique => $lane, lane => $lane });
        my %return = (lane => $vrlane);
        
        my ($vrlib) = $vrlane->related(incoming => { namespace => 'VRTrack', label => 'Library' });
        if ($vrlib) {
            if ($vrlib->id ne $library) {
                if ($enforce) {
                    $vrlib = $self->add('Library', { id => $library });
                    $vrlib->relate_to($vrlane, 'sequenced', selfish => 1);
                }
                else {
                    return \%return;
                }
            }
        }
        else {
            $vrlib = $self->add('Library', { id => $library }, outgoing => { type => 'sequenced', node => $vrlane });
        }
        $return{library} = $vrlib;
        
        my ($vrsam) = $vrlib->related(incoming => { namespace => 'VRTrack', label => 'Sample' });
        if ($vrsam) {
            if ($vrsam->name ne $sample) {
                if ($enforce) {
                    $vrsam = $self->add('Sample', { name => $sample });
                    $vrsam->relate_to($vrlib, 'prepared', selfish => 1);
                }
                else {
                    return \%return;
                }
            }
        }
        else {
            $vrsam = $self->add('Sample', { name => $sample }, outgoing => { type => 'prepared', node => $vrlib });
        }
        $return{sample} = $vrsam;
        
        # because a sample can be associated with more than 1 study, and study
        # is optional, we just set it if provided and don't return existing
        # studies otherwise
        if ($study) {
            $return{study} = $self->add('Study', { id => $study }, outgoing => { type => 'member', node => $vrsam });
            
            # and since a study can be in more than 1 group, we ignore group
            # if they didn't supply the study to attach it to
            if ($group) {
                $return{group} = $self->add('Group', { name => $group }, outgoing => { type => 'has', node => $return{study} });
            }
        }
        
        if ($taxon) {
            my ($vrtax) = $vrsam->related(incoming => { namespace => 'VRTrack', label => 'Taxon' });
            if ($vrtax) {
                if ($vrtax->id ne $taxon) {
                    if ($enforce) {
                        $vrtax = $self->add('Taxon', { id => $taxon });
                        $vrtax->relate_to($vrsam, 'member', selfish => 1);
                    }
                    else {
                        return \%return;
                    }
                }
            }
            else {
                $vrtax = $self->add('Taxon', { id => $taxon }, outgoing => { type => 'member', node => $vrsam });
            }
            $return{taxon} = $vrtax;
        }
        
        return \%return;
    }
    
    method add_file (Str $path) {
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        return $vrpipe_schema->path_to_filesystemelement($path);
    }
    
    method get_file (Str $path) {
        $vrpipe_schema ||= VRPipe::Schema->create('VRPipe');
        return $vrpipe_schema->path_to_filesystemelement($path, only_get => 1);
    }
    
    # donor's just have meaningless ids; user will want to see the
    # donor's control sample public name and the most recent created_date of its
    # member samples, so we'll need some cypher that get the donor's samples
    # and the relationship between them. You're supposed to append the return
    # value of this to some cypher that first matches your desired donor(s).
    method donor_to_sample_match_cypher (Str $donor_identifer = 'donor') {
        my $sample_labels = $self->cypher_labels('Sample');
        return "MATCH ($donor_identifer)-[ds_rel]->(d_sample:$sample_labels) RETURN $donor_identifer,ds_rel,d_sample";
    }
    
    # if user ran a cypher query based on donor_to_sample_match_cypher()
    # we have sample nodes related to our donor nodes; figure out
    # which samples belong to which donors, then add the relevant
    # sample properties to the donor node. Input is the output of run_cypher().
    method add_sample_info_to_donors (HashRef $graph_data) {
        my %rels;
        foreach my $rel (@{ $graph_data->{relationships} }) {
            push(@{ $rels{ $rel->{startNode} } }, $rel->{endNode});
        }
        
        my %nodes;
        foreach my $node (@{ $graph_data->{nodes} }) {
            $nodes{ $node->{label} }->{ $node->{id} } = $node;
        }
        
        foreach my $donor (values $nodes{Donor}) {
            my ($most_recent_date, @controls, $shortest);
            foreach my $sample_id (@{ $rels{ $donor->{id} } || next }) {
                my $s_props = $nodes{Sample}->{$sample_id}->{properties};
                
                my $cd = $s_props->{created_date};
                if ($cd && (!defined $most_recent_date || $cd > $most_recent_date)) {
                    $most_recent_date = $cd;
                }
                
                my $name = $s_props->{public_name} || $s_props->{name};
                if (!defined $shortest || length($name) < length($shortest)) {
                    $shortest = $name;
                }
                
                push(@controls, $name) if $s_props->{control};
            }
            
            if (@controls || $shortest) {
                $donor->{properties}->{example_sample} = @controls == 1 ? $controls[0] : $shortest;
            }
            if ($most_recent_date) {
                # convert to yyyy-mm-dd
                my $dt = DateTime->from_epoch(epoch => $most_recent_date);
                $donor->{properties}->{last_sample_added_date} = $dt->ymd;
            }
        }
        
        return [sort { $b->{properties}->{last_sample_added_date} cmp $a->{properties}->{last_sample_added_date} || $a->{properties}->{example_sample} cmp $b->{properties}->{example_sample} } values $nodes{Donor}];
    }
}

1;
