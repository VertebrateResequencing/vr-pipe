
=head1 NAME

VRPipe::Steps::bam_metadata_many - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::bam_metadata_many extends VRPipe::Steps::bam_metadata {
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            my $bamcheck_exe = $options->{bamcheck_exe};
            my $store_pg     = $options->{store_original_pg_chain};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my (@outfiles, @io_map, @pg_inputs);
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $ifile = $bam_file->path;
                
                # run bamcheck if we don't have enough metadata
                my $meta       = $bam_file->metadata;
                my $meta_count = 0;
                foreach my $type (@{ $self->meta_to_check }) {
                    $meta_count++ if $meta->{$type};
                }
                unless ($meta_count == @{ $self->meta_to_check }) {
                    my $check_file = $self->output_file(basename => $ifile->basename . '.bamcheck', type => 'txt', temporary => 1);
                    push @io_map,   $bam_file->id . '=>' . $check_file->id;
                    push @outfiles, $check_file;
                }
                
                # we'll also check the header for existing PG lines and store
                # those as metadata
                if ($store_pg && !defined $meta->{original_pg_chain}) {
                    push @pg_inputs, $bam_file->id;
                }
            }
            my $files_string = join ',', @io_map;
            my $this_cmd = "use VRPipe::Steps::bam_metadata_many; VRPipe::Steps::bam_metadata_many->stats_from_bamcheck_many(bamcheck => $bamcheck_exe, files => { $files_string });";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
            
            if (@pg_inputs) {
                $self->dispatch_vrpipecode("use VRPipe::Steps::bam_metadata_many; VRPipe::Steps::bam_metadata_many->store_pg_chain_many(bams => [qw(@pg_inputs)]);", $req);
            }
        };
    }
    
    method description {
        return "Takes a bam file and associates metadata with the file in the VRPipe database, making the bam file usable in other bam-related Steps";
    }
    
    method stats_from_bamcheck_many (ClassName|Object $self: Str :$bamcheck, HashRef :$files) {
        keys %{$files} > 0 || $self->throw("You must supply file ids to this method");
        while (my ($input, $output) = each %{$files}) {
            my $input_path  = VRPipe::File->get(id => $input)->path;
            my $output_path = VRPipe::File->get(id => $output)->path;
            my $cmd_line    = qq[$bamcheck $input_path > $output_path];
            my $ok          = $self->stats_from_bamcheck($cmd_line);
            unless ($ok) {
                $self->throw("cmd [$cmd_line] failed");
            }
        }
        return 1;
    }
    
    method store_pg_chain_many (ClassName|Object $self: ArrayRef :$bams) {
        scalar @{$bams} > 0 || $self->throw("You must supply file ids to this method");
        foreach my $id (@{$bams}) {
            my $input_path = VRPipe::File->get(id => $id)->path;
            my $ok = $self->store_pg_chain(bam => qq[$input_path]);
            unless ($ok) {
                $self->throw("store_pg_chain failed for $input_path");
            }
        }
        return 1;
    }
}

1;
