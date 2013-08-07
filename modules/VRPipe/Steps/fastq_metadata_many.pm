
=head1 NAME

VRPipe::Steps::fastq_metadata_many - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::fastq_metadata_many extends VRPipe::Steps::fastq_metadata {
    method body_sub {
        return sub {
            my $self = shift;
            
            my $fastqcheck_exe = $self->options->{fastqcheck_exe};
            
            my (@outfiles, @io_map);
            foreach my $fq_file (@{ $self->inputs->{fastq_files} }) {
                my $ifile = $fq_file->path;
                
                # we expect the datasource to fill in most if not all metadata,
                # but some things may need calculating
                my $meta = $fq_file->metadata;
                unless ($meta->{bases} && $meta->{reads} && $meta->{avg_read_length}) {
                    # run fastqcheck to generate stats
                    # record the fact that we (probably) made the fastqcheck
                    # file, even though it isn't in our outputs_definition *** or should we be able to have optional outputs?
                    my $fqc_file = $self->output_file(basename => $ifile->basename . '.fastqcheck', type => 'txt', temporary => 1);
                    push @io_map,   $fq_file->id . '=>' . $fqc_file->id;
                    push @outfiles, $fqc_file;
                }
            }
            if (@io_map) {
                my $files_string = join ',', @io_map;
                my $req = $self->new_requirements(memory => 500, time => 1);
                my $this_cmd = "use VRPipe::Steps::fastq_metadata_many; VRPipe::Steps::fastq_metadata_many->stats_from_fastqcheck_many(fastqcheck => $fastqcheck_exe, files => { $files_string });";
                $self->dispatch_vrpipecode($this_cmd, $req);
            }
        };
    }
    
    method description {
        return "Takes a fastq file and associates metadata with the file in the VRPipe database, making the fastq file usable in other fastq-related Steps";
    }
    
    method post_process_sub {
        return sub {
            my $self = shift;
            
            my $all_ok = 1;
            foreach my $ofile (@{ $self->inputs->{fastq_files} }) {
                my $meta = $ofile->metadata;
                unless ($meta->{bases} && $meta->{reads} && $meta->{avg_read_length}) {
                    $self->warn("some metadata was missing from the dataelement of " . $ofile->path . ", and our fastqcheck run failed to resolved that");
                    $all_ok = 0;
                }
            }
            
            return $all_ok;
        };
    }
    
    method stats_from_fastqcheck_many (ClassName|Object $self: Str :$fastqcheck, HashRef :$files) {
        keys %{$files} > 0 || $self->throw("You must supply file ids to this method");
        while (my ($input, $output) = each %{$files}) {
            my $input_path  = VRPipe::File->get(id => $input)->path;
            my $output_path = VRPipe::File->get(id => $output)->path;
            my $cmd_line    = qq[$fastqcheck $input_path > $output_path];
            my $ok          = $self->stats_from_fastqcheck($cmd_line);
            unless ($ok) {
                $self->throw("cmd [$cmd_line] failed");
            }
        }
        return 1;
    }
}

1;
