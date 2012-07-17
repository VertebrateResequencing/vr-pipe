=head1 NAME

VRPipe::Steps::md5_file_production - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::md5_file_production with VRPipe::StepRole {
    method options_definition {
        return { md5_files_in_source_dir => VRPipe::StepOption->create(description => 'should .md5 files be created in the same directory as their parent? (otherwise they go to the pipeline output directory)',
                                                                    optional => 1,
                                                                    default_value => 0) };
    }
    method inputs_definition {
        return { md5_file_input => VRPipe::StepIODefinition->create(type => 'any', max_files => -1, description => 'one or more files you want to calculate the md5 checksum of') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $input_files = $self->inputs->{md5_file_input};
            my $md5_files_in_source_dir = $self->options->{md5_files_in_source_dir};
            
            @$input_files || return 1;
            foreach my $vrfile (@$input_files) {
               my $ifile = $vrfile->path;
               my $ofile = $self->output_file(output_key => 'md5_files',
                                              $md5_files_in_source_dir ? (output_dir => $ifile->dir) : (),
                                              basename => $ifile->basename.'.md5',
                                              type => 'txt',
                                              metadata => {md5_checksum_of => $ifile->stringify})->path;
               $self->dispatch([qq{md5sum $ifile > $ofile}, $self->new_requirements(memory => 500, time => 1)]);
            }
        };
    }
    method outputs_definition {
        return { md5_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => '.md5 files for each input file containing its md5 checksum',
                                                            metadata => {md5_checksum_of => 'absolute path of the file this .md5 file was generated for'}) };
    }
    method post_process_sub {
        return sub {
            my $self = shift;
            my $output_files = $self->outputs->{md5_files};
            foreach my $ofile (@$output_files) {
               my $content = $ofile->slurp;
               $content || return 0;
               my ($md5) = split(" ", $content);
               my $ifile = VRPipe::File->get(path => $ofile->metadata->{md5_checksum_of});
               $ifile->md5($md5);
               $ifile->update;
            }
            return 1;
        };
    }
    method description {
        return "Takes a file, calculates its md5 checksum, produces a file called <input filename>.md5, and updates the persistent database with the md5 of the file";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
