
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
        return { iget => 'iget', ichksum => 'ichksum', imeta => 'imeta' };
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
            my $imeta      = $opts->{imeta_exe};
            my $dir        = $opts->{vrlane_storage_dir};
            my $get_inputs = $opts->{irods_download_input_files};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $file (@{ $self->inputs->{files} }) {
                my $source_file = $file->path->stringify;
                my $source_meta = $file->metadata;
                delete $source_meta->{expected_md5};
                delete $source_meta->{md5};
                
                my $afiles = delete $source_meta->{irods_analysis_files};
                my @afiles = ref($afiles) ? @$afiles : ($afiles) if $afiles;
                my $irods_local_storage_dir;
                foreach my $afile (@afiles) {
                    $irods_local_storage_dir = $dir;
                    my $download_path = file($dir, $afile);
                    my $download_file = VRPipe::File->create(path => $download_path);
                    
                    $self->output_file(output_key => 'analysis_files', output_dir => $download_file->dir, basename => $download_file->basename, type => 'any');
                    my $md5 = $download_file->meta_value('md5');
                    $download_file->merge_metadata({ source_file => $source_file, %{$source_meta} });
                    $download_file->add_metadata({ irods_path => $afile, $md5 ? (md5 => $md5) : () }, replace_data => 1);
                    
                    next if $download_file->s;
                    
                    $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_analysis_files_download; VRPipe::Steps::irods_analysis_files_download->get_file(source => q[$afile], dest => q[$download_path], iget => q[$iget], ichksum => q[$ichksum], imeta => q[$imeta], add_metadata => 1);], $req, { output_files => [$download_file], block_and_skip_if_ok => 1 });
                }
                
                if ($get_inputs) {
                    my $meta = $file->metadata;
                    $meta->{irods_local_storage_dir} = $irods_local_storage_dir if $irods_local_storage_dir;
                    
                    if (!$file->s) {
                        # actually download to the path of input ($file)
                        my $dest = $self->output_file(output_key => 'input_files', output_dir => $file->dir, basename => $file->basename, type => $file->type, metadata => $meta)->path;
                        $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_analysis_files_download; VRPipe::Steps::irods_analysis_files_download->get_file(source => q[$meta->{irods_path}], dest => q[$dest], iget => q[$iget], ichksum => q[$ichksum]);], $req);
                    }
                    else {
                        # check it wasn't us that previously downloaded $file
                        my %output_by = map { $_->id => 1 } $file->output_by;
                        unless (exists $output_by{ $self->step_state->id }) {
                            # symlink our existing input file to the pipeline
                            # output dir so that if this step is restarted, we
                            # won't delete our input file
                            my $ofile = $self->output_file(output_key => 'input_files', basename => $file->basename, type => $file->type, metadata => $file->metadata);
                            $file->symlink($ofile);
                        }
                    }
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            analysis_files => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'a file on a local disc, extracted from iRODs',
                min_files   => 0,
                max_files   => -1
            ),
            input_files => VRPipe::StepIODefinition->create(
                type            => 'any',
                description     => 'a file on a local disc, extracted from iRODs',
                min_files       => 0,
                max_files       => -1,
                check_existence => 0
            )
        };
    }
    
    method description {
        return "Get files out of iRODs and store them on local disc. The file paths to download are found in the input file irods_analysis_files metadata. Optionally also download the input file.";
    }
}

1;
