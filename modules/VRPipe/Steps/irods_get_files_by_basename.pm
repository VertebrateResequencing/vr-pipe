
=head1 NAME

VRPipe::Steps::irods_get_files_by_basename - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012,2014,2015 Genome Research Limited.

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
    use VRPipe::Schema;
    
    our $schema;
    
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
            ),
            irods_convert_cram_to_bam => VRPipe::StepOption->create(
                description   => "when getting a cram file from irods, convert it to a bam file; 0 turns this off, the absolute path to a samtools v1+ executable turns this on",
                optional      => 1,
                default_value => 0
            ),
            iget_args => VRPipe::StepOption->create(
                description   => 'optional argument to iget when retrieving files',
                optional      => 1,
                default_value => '-K -f'
            ),
            ichksum_args => VRPipe::StepOption->create(
                description => 'optional argument to ichksum when retrieving files',
                optional    => 1,
            ),
            local_root_dir => VRPipe::StepOption->create(
                description => 'optional root directory to store the irods files; required if you did not set this option in your irods datasource, ignored otherwise',
                optional    => 1,
            )
        };
    }
    
    method inputs_definition {
        return {
            basenames => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'either a) file paths that do not exist yet - the basename will be used to find the file in iRODs, and the file will be saved at this full path, or b) irods absolute file paths of files that exist in irods (set with the irods: protocol)',
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
            
            my $opts             = $self->options;
            my $zone             = $opts->{irods_get_zone};
            my $iget             = $opts->{iget_exe};
            my $iget_args        = $opts->{iget_args} || '';
            my $ichksum_args     = $opts->{ichksum_args} || '';
            my $iquest           = $opts->{iquest_exe};
            my $ichksum          = $opts->{ichksum_exe};
            my $samtools         = $opts->{irods_convert_cram_to_bam};
            my $cram_to_bam_mode = $samtools && -x $samtools;
            my $local_root_dir   = $opts->{local_root_dir};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $file (@{ $self->inputs->{basenames} }) {
                my $basename        = $file->basename;
                my $meta            = $file->metadata;
                my $output_basename = $basename;
                my $output_type     = $file->type;
                my $extra           = '';
                if ($cram_to_bam_mode && $basename =~ /\.cram$/) {
                    $output_basename =~ s/\.cram$/.bam/;
                    $output_type = 'bam';
                    $extra       = ", samtools_for_cram_to_bam => q[$samtools]";
                }
                
                my ($output_dir, $irods_abs_path);
                if ($file->pro eq 'irods') {
                    # we'll download to the path of the irods file with
                    # local_root_dir as our root dir
                    $local_root_dir || $self->throw("local_root_dir must be specified for irods:/ files");
                    $irods_abs_path = $file->protocolless_path;
                    $output_dir = dir($local_root_dir, $irods_abs_path->dir);
                }
                else {
                    # we'll download to the path of our input file
                    $output_dir     = $file->dir;
                    $irods_abs_path = $meta->{irods_path};
                }
                
                my $out_abs_path = file($output_dir, $output_basename);
                
                if (!-s $out_abs_path) {
                    my $dest = $self->output_file(output_key => 'local_files', output_dir => $output_dir, basename => $output_basename, type => $output_type, metadata => $meta);
                    
                    if ($iget_args) {
                        $extra .= ", iget_args => q[$iget_args]";
                    }
                    if ($ichksum_args) {
                        $extra .= ", ichksum_args => q[$ichksum_args]";
                    }
                    
                    # if we have the full irods path, get the file directly,
                    # otherwise search for it by basename
                    if ($irods_abs_path) {
                        $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_get_files_by_basename; VRPipe::Steps::irods_get_files_by_basename->get_file(source => q[$irods_abs_path], dest => q[$out_abs_path], iget => q[$iget], ichksum => q[$ichksum]$extra);], $req);
                    }
                    else {
                        $self->dispatch_vrpipecode(qq[use VRPipe::Steps::irods_get_files_by_basename; VRPipe::Steps::irods_get_files_by_basename->get_file_by_basename(basename => q[$basename], dest => q[$out_abs_path], zone => q[$zone], iget => q[$iget], iquest => q[$iquest], ichksum => q[$ichksum]$extra);], $req);
                    }
                }
                else {
                    # symlink our existing output file to the pipeline output
                    # dir so that if this step is restarted, we won't delete it
                    my $ofile = $self->output_file(output_key => 'local_files', basename => $output_basename, type => $output_type, metadata => $meta);
                    
                    # relate source file to dest file in the graph database
                    $schema ||= VRPipe::Schema->create('VRPipe');
                    my $source_graph_node = $schema->get('File', { path => $file->protocolless_path, protocol => $file->protocol });
                    if ($source_graph_node) {
                        $self->relate_input_to_output($source_graph_node, 'imported', $out_abs_path->stringify);
                    }
                    
                    VRPipe::File->create(path => $out_abs_path)->symlink($ofile);
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
        return "Get files out of iRODs and store them on local disc. Tries to find the files in iRODs based on the basenames of the input files, if you do not supply the full iRODs path.";
    }
}

1;
