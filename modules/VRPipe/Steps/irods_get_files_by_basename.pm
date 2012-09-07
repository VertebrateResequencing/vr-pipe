
=head1 NAME

VRPipe::Steps::irods_get_files_by_basename - a step

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

class VRPipe::Steps::irods_get_files_by_basename extends VRPipe::Steps::irods {
    method _build_irods_exes {
        return { iget => 'iget', iquest => 'iquest', ichksum => 'ichksum' };
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            irods_get_zone => VRPipe::StepOption->create(
                description   => 'the zone (top level directory) where your data is stored in iRODs',
                optional      => 1,
                default_value => 'seq'
            )
        };
    }
    
    method inputs_definition {
        return {
            basenames => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'file paths that do not exist yet - the basename will be used to find the file in iRODs, and the file will be saved at this full path',
                max_files   => -1,
                metadata    => {
                    expected_md5 => 'the md5 checksum the file is supposed to have',
                    optional     => ['expected_md5']
                },
                check_existence => 0
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $opts    = $self->options;
            my $zone    = $opts->{irods_get_zone};
            my $iget    = $opts->{iget_exe};
            my $iquest  = $opts->{iquest_exe};
            my $ichksum = $opts->{ichksum_exe};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $file (@{ $self->inputs->{basenames} }) {
                if (!$file->s) {
                    # download to the path of our input file, which doesn't
                    # exist yet
                    my $basename = $file->basename;
                    my $dest = $self->output_file(output_key => 'local_files', output_dir => $file->dir, basename => $basename, type => $file->type, metadata => $file->metadata)->path;
                    $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_get_files_by_basename; VRPipe::Steps::irods_get_files_by_basename->get_file_by_basename(basename => q[$basename], dest => q[$dest], zone => q[$zone], iget => q[$iget], iquest => q[$iquest], ichksum => q[$ichksum]);], $req);
                }
                else {
                    # symlink our existing input file to the pipeline output dir
                    # so that if this step is restarted, we won't delete our
                    # input file
                    my $ofile = $self->output_file(output_key => 'local_files', basename => $file->basename, type => $file->type, metadata => $file->metadata);
                    $file->symlink($ofile);
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            local_files => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'a file on a local disc, extracted from iRODs',
                max_files   => -1
            )
        };
    }
    
    method description {
        return "Get files out of iRODs and store them on local disc. Tries to find the files in iRODs based on the basenames of the input files.";
    }
}

1;
