
=head1 NAME

VRPipe::Steps::archive_index_files - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::archive_index_files with VRPipe::StepRole {
    use Digest::MD5;
    
    method options_definition {
        return {};
    }
    
    method inputs_definition {
        return {
            file  => VRPipe::StepIODefinition->create(type => 'any', description => 'a file that has already been archived'),
            index => VRPipe::StepIODefinition->create(type => 'any', description => 'an index file associated with file that should be archived to the same directory')
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            # we do not allow for multiple input files because if 1 gets moved
            # successfully but 1 does not and the user does a reset to fix
            # the failure, the successfully moved file gets deleted and our
            # input no longer exists!
            my ($file)  = @{ $self->inputs->{file} };
            my ($index) = @{ $self->inputs->{index} };
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my $ofile = VRPipe::File->create(
                path     => file($file->resolve->dir, $index->basename),
                type     => $index->type,
                metadata => $index->metadata
            )->path;
            
            # VRPipe::File->move will do a copy, check md5s, and only then
            # delete the source, so it is nice and safe
            my $ifile = $index->path;
            $self->dispatch_vrpipecode("VRPipe::File->get(path => q[$ifile])->move(VRPipe::File->get(path => q[$ofile]), check_md5s => 1) || die q[move failed];", $req);
        };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Move an index file to the same archive directory as it's associated file";
    }
    
    method max_simultaneous {
        return 0;            # archive files are small, so don't bother with a limit on these
    }
}

1;
