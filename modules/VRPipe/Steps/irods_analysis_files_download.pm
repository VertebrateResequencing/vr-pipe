
=head1 NAME

VRPipe::Steps::irods_analysis_files_download - a step

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

class VRPipe::Steps::irods_analysis_files_download extends VRPipe::Steps::irods {
    method _build_irods_exes {
        return { iget => 'iget', ichksum => 'ichksum' };
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            vrlane_storage_dir         => VRPipe::StepOption->create(description => 'the absolute path to the base directory that lane-related files will be stored in (unnecessary for sequencing data)', optional => 1),
            irods_download_input_files => VRPipe::StepOption->create(description => 'boolean; optionally also download the input file',                                                                    optional => 1, default_value => 0)
        };
    }
    
    method inputs_definition {
        return {
            files => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'files that have irods_analysis_files metadata',
                max_files   => -1,
                metadata    => {
                    irods_analysis_files => '1 or more paths in irods to download',
                    irods_path           => 'path to the input file in irods',
                    optional             => ['irods_analysis_files']               # it's optional to let people use a pipeline containing this step when they have no irods_analysis_files
                },
                check_existence => 0
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $opts       = $self->options;
            my $iget       = $opts->{iget_exe};
            my $ichksum    = $opts->{ichksum_exe};
            my $dir        = $opts->{vrlane_storage_dir};
            my $get_inputs = $opts->{irods_download_input_files};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $file (@{ $self->inputs->{files} }) {
                my @afiles = $file->meta_value("irods_analysis_files");
                foreach my $afile (@afiles) {
                    my $download_path = file($dir, $afile);
                    my $download_file = VRPipe::File->create(path => $download_path);
                    next if $download_file->s;
                    
                    $self->output_file(output_key => 'local_files', output_dir => $download_file->dir, basename => $download_file->basename, type => 'any', metadata => { source_file => $file->path->stringify })->path;
                    
                    $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_analysis_files_download; VRPipe::Steps::irods_analysis_files_download->get_file(source => q[$afile], dest => q[$download_path], iget => q[$iget], ichksum => q[$ichksum]);], $req, { output_files => [$download_file], block_and_skip_if_ok => 1 });
                }
                
                if ($get_inputs) {
                    # (pretty much the same as the irods_get_files_by_basename
                    #  step)
                    if (!$file->s) {
                        # download to the path of our input file, which doesn't
                        # exist yet
                        my $meta = $file->metadata;
                        my $dest = $self->output_file(output_key => 'local_files', output_dir => $file->dir, basename => $file->basename, type => $file->type, metadata => $meta)->path;
                        
                        $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_analysis_files_download; VRPipe::Steps::irods_analysis_files_download->get_file(source => q[$meta->{irods_path}], dest => q[$dest], iget => q[$iget], ichksum => q[$ichksum]);], $req);
                    }
                    else {
                        # symlink our existing input file to the pipeline output dir
                        # so that if this step is restarted, we won't delete our
                        # input file
                        my $ofile = $self->output_file(output_key => 'local_files', basename => $file->basename, type => $file->type, metadata => $file->metadata);
                        $file->symlink($ofile);
                    }
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            local_files => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'a file on a local disc, extracted from iRODs',
                min_files   => 0,
                max_files   => -1
            )
        };
    }
    
    method description {
        return "Get files out of iRODs and store them on local disc. The file paths to download are found in the input file irods_analysis_files metadata. Optionally also download the input file.";
    }
}

1;
