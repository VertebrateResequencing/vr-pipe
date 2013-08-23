
=head1 NAME

VRPipe::Steps::vrtrack_update_mapstats_hipsci - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Steps::vrtrack_update_mapstats_hipsci extends VRPipe::Steps::vrtrack_update {
    around options_definition {
        return { %{ $self->$orig } };
    }
    
    method inputs_definition {
        return {
            diff_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Bed file containing the diff of control and stem cell bed files',
                max_files   => -1,                                                               # -1 = As many as you like
                min_files   => 0
            )
        };
    }
    
    method body_sub {
        return sub {
            #add metadata for analysis type to add to mapstats, i.e. penncnv/quantisnp.....
            my $self = shift;
            my $opts = $self->options;
            my $db   = $opts->{vrtrack_db};
            my $req  = $self->new_requirements(memory => 500, time => 1);
            foreach my $diff_file (@{ $self->inputs->{diff_files} }) {
                my $diff_path = $diff_file->path;
                my $lane      = $diff_file->metadata->{lane};
                if (defined($diff_file->metadata->{cnv_total})) {
                    my $cmd = "use VRPipe::Steps::vrtrack_update_mapstats_hipsci; VRPipe::Steps::vrtrack_update_mapstats_hipsci->update_mapstats(db => q[$db], diff => q[$diff_path], lane => q[$lane]);";
                    $self->dispatch_vrpipecode($cmd, $req);
                }
            }
        };
    }
    
    method description {
        return "Add the CNV intersection and diff statistics to the mapstatst for hipsci sampless to the VRTrack database, so that they're accessible with QCGrind etc.";
    }
    
    method update_mapstats (ClassName|Object $self: Str :$db!, Str|File :$diff!, Str :$lane!) {
        my $diff_file = VRPipe::File->get(path => $diff);
        my $meta = $diff_file->metadata;
        $diff_file->disconnect;
        my $cnv_total         = $meta->{cnv_total};
        my $cnv_diff          = $meta->{cnv_minus_control};
        my $cnv_analysis_type = $meta->{cnv_analysis_type};
        
        # get the lane and mapstats object from VRTrack
        my $vrtrack = $self->get_vrtrack(db => $db);
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("No lane named '$lane' in database '$db'");
        
        $diff_file->update_stats_from_disc;
        my $fh = $diff_file->openr;
        while (my $coordinates = <$fh>) {
            chomp $coordinates;
            $coordinates =~ s/\t/:/g;
            
            my $worked = $vrtrack->transaction(
                sub {
                    my $mapstats = $vrlane->add_mapping();
                    $mapstats->raw_bases($cnv_total);
                    $mapstats->clip_bases($cnv_diff);
                    $mapstats->genotype_found($coordinates);
                    $mapstats->genotype_expected($cnv_analysis_type);
                    $mapstats->update;
                }
            );
            
            unless ($worked) {
                $self->throw($vrtrack->{transaction_error});
            }
        }
    }
}

1;
