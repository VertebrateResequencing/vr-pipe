=head1 NAME

VRPipe::Steps::vrtrack_update_improved - a step

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

class VRPipe::Steps::vrtrack_update_improved extends VRPipe::Steps::vrtrack_update {
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', 
                                                            description => 'bam file that has been improved', 
                                                            max_files => -1,
                                                            metadata => {lane => 'lane name (a unique identifer for this sequencing run, aka read group)'}) };
    }
    method body_sub {
        return sub { 
            my $self = shift;
	    my $db = $self->options->{vrtrack_db};
            my $req = $self->new_requirements(memory => 500, time => 1);
	    
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam_file->path;
		my $lane = $bam_file->metadata->{lane};
                my $cmd = "use VRPipe::Steps::vrtrack_update_improved; VRPipe::Steps::vrtrack_update_improved->update_improved(db => q[$db], bam => q[$bam_path], lane => q[$lane]);";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    method description {
        return "Add the improved bam to the VRTrack file table for the lane.";
    }

    method update_improved (ClassName|Object $self: Str :$db!, Str|File :$bam!, Str :$lane!) {
	my $bam_file = VRPipe::File->get(path => $bam);  	 
	my $meta = $bam_file->metadata;
	$bam_file->disconnect;
	
	# get the lane object from VRTrack
	my $vrtrack = $self->get_vrtrack(db => $db);
	my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("No lane named '$lane' in database '$db'");
	
	my $worked = $vrtrack->transaction(sub {
	    # get/create a new file entry for us
	    my $vrfile_name = 'VRPipe::File::'.$bam_file->id;
	    my $vrfile = $vrlane->get_file_by_name($vrfile_name); # this shouldn't really ever work
	    unless ($vrfile) {
		$vrfile = $vrlane->add_file($vrfile_name);
	    }
	    $vrfile->type(5);
	    $vrfile->is_processed(import => 1);
	    $vrfile->md5($bam_file->md5);
	    $vrfile->update;
	    $vrfile->is_processed(mapped => 1);
	    $vrfile->update;
	    $vrfile->is_processed(improved => 1);
	    $vrfile->update;
	    
	    # also update the lane; in case previously we didn't fill out reads and
	    # bases (eg. we ran bam import and QC in exome intervals only), we'll do
	    # that now as well
	    $vrlane->raw_bases($meta->{bases}) if $meta->{bases};
	    $vrlane->raw_reads($meta->{reads}) if $meta->{reads};
	    $vrlane->is_paired($meta->{paired} ? 1 : 0) if exists $meta->{paired};
	    $vrlane->read_len(int($meta->{avg_read_length})) if $meta->{avg_read_length};
	    $vrlane->is_processed(mapped => 1);
	    $vrlane->update;
	    $vrlane->is_processed(improved => 1);
	    $vrlane->update;
	});
	
	unless ($worked) {
	    $self->throw($vrtrack->{transaction_error});
	}
    }
}

1;
