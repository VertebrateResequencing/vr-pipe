
=head1 NAME

VRPipe::Steps::archive_files - a step

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

class VRPipe::Steps::archive_files with VRPipe::StepRole {
    use Digest::MD5;
    
    method options_definition {
        return { disc_pool_file => VRPipe::StepOption->create(description => 'path to a file with an absolute path of a storage root directory on each line') };
    }
    
    method inputs_definition {
        return { file => VRPipe::StepIODefinition->create(type => 'any', description => 'a file that should be archived') };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options        = $self->options;
            my $disc_pool_file = file($options->{disc_pool_file});
            $self->throw("disc_pool_file must be an absolute path") unless $disc_pool_file->is_absolute;
            
            # parse the disc_pool_file; we don't cache these results in a class
            # variable because disc_pool_file could be updated at any time.
            my $fh = $disc_pool_file->openr;
            my @dirs;
            while (<$fh>) {
                chomp;
                next if /^#/;
                my $dir = $_;
                next unless dir($dir)->is_absolute;
                push(@dirs, $dir);
            }
            my $max = @dirs;
            
            my $req = $self->new_requirements(memory => 500, time => 1); # md5sum can be memory intensive??
            
            # we do not allow for multiple input files because if 1 gets moved
            # successfully but 1 does not and the user does a reset to fix
            # the failure, the successfully moved file gets deleted and our
            # input no longer exists!
            my ($file) = @{ $self->inputs->{file} };
            # pick a random directory from the pool
            my $root_dir = $dirs[int(rand($max))];
            
            # pick subdirs based on md5sum of $file's id
            my $dmd5     = Digest::MD5->new();
            my $ifile_id = $file->id;
            $dmd5->add($ifile_id);
            my $md5 = $dmd5->hexdigest;
            my @chars = split("", $md5);
            
            my $ofile = $self->output_file(
                output_key => 'moved_file',
                basename   => $file->basename,
                output_dir => dir($root_dir, @chars[0 .. 3], $ifile_id),
                type       => $file->type,
                metadata   => $file->metadata
            )->path;
            
            # VRPipe::File->move will do a copy, check md5s, and only then
            # delete the source, so it is nice and safe
            my $ifile = $file->path;
            $self->dispatch_vrpipecode("VRPipe::File->get(path => q[$ifile])->move(VRPipe::File->get(path => q[$ofile]), check_md5s => 1);", $req);
        };
    }
    
    method outputs_definition {
        return {
            moved_file => VRPipe::StepIODefinition->create(
                type        => 'any',
                description => 'the input file at its new location'
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Safely move a file from one disc to a pool of one or more other discs, eg. for archival purposes";
    }
    
    method max_simultaneous {
        return 25;
    }
}

1;
