
=head1 NAME

VRPipe::Steps::vcf_to_irods - a step

=head1 DESCRIPTION

Adds vcf files and associated metadata to irods. Requires that the vcfs have sample(s) and md5 checksum vrpipe metadata.
The vcfs are added under an irods root collection using a directory structure constructed from the md5 checksum.
Also optionally adds the vcf index to irods.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::vcf_to_irods with VRPipe::StepRole {
    
    method options_definition {
        return {
            'irods_root' => VRPipe::StepOption->create( description  => 'irods root path to which to add objects'),
            'has_index' => VRPipe::StepOption->create( description  => 'boolean, indicates vcf has an index which needs also adding to irods', optional => 1, default_value => 1),
            'study' => VRPipe::StepOption->create( description  => 'study to which vcf belongs, to be added as irods metadata, if not found in vcf metadata', optional => 1),
            'release' => VRPipe::StepOption->create( description  => 'release to which vcf belongs, to be added irods as metadata, if not found in vcf metadata', optional => 1),
        };
    }

    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'vcf files',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options = $self->options;
            my $irods_root = $options->{irods_root};
            my $has_index = $options->{has_index};
            my $study = $options->{study};
            my $release = $options->{release};
            
            my $req = $self->new_requirements(memory => 100, time => 1);

            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {

                my $input_path = $vcf_file->path;
                my $basename = $vcf_file->basename;

                my $cmd = "use VRPipe::Steps::vcf_to_irods; VRPipe::Steps::vcf_to_irods->add_vcf_to_irods('$input_path', '$irods_root', '$has_index','$study', '$release');";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    
    method outputs_definition {
        return { };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Takes a vcf file and adds the file and its assocated metadata to irods";
    }
    
    method add_vcf_to_irods (ClassName|Object $self: Str|File $input_vcf_path, Str $irods_root, Bool $has_index , Str $study, Str $release) {

        my $vcf_file = VRPipe::File->get(path => $input_vcf_path);

        # Add collection using directory structure constructed from md5 checksum, eg
        # eg 198f9d136fcd9b7f83b3bf32dc25f181: /uk10k/uk10k/STUDY/1/9/8/f/d136fcd9b7f83b3bf32dc25f181

        my $vcf_meta = $vcf_file->metadata;
        my $md5 = $vcf_meta->{md5};
        $self->throw("Could not get md5 checksum metadata for $input_vcf_path") unless $md5;

        my ($md51,$md52,$md53,$md54) = split(//,$md5);
        my $md55=substr($md5,5);
        my $idir="$irods_root/$md51/$md52/$md53/$md54/$md55";

        #add collection and object
        system("imkdir", "-p", "$idir") && $self->throw("failed to imkdir $idir");
        system("iput", "$input_vcf_path", "$idir") && $self->throw("failed to iput $input_vcf_path to $idir");
        if ($has_index) {
            system("iput", "$input_vcf_path.tbi", "$idir") && $self->throw("failed to iput $input_vcf_path.tbi to $idir");
        }

        my $basename = $vcf_file->basename;
        system("ils", "-r", "$idir/$basename") && $self->throw("failed to ils $idir/$basename");

        # Add metadata
        my @cmd = ("imeta", "add", "-d", "$idir/$basename", "type", "vcf");
        system(@cmd) && $self->throw("failed to run '@cmd'");

        @cmd = ("imeta", "add", "-d", "$idir/$basename", "md5", "$md5");
        system(@cmd) && $self->throw("failed to run '@cmd'");

        my $meta_count=2;

        # sample names are catanated into a string and delimited by '#'
        my @samples = split(/#/,$vcf_meta->{samples});
        foreach my $sample(@samples) {
            @cmd = ("imeta", "add", "-d", "$idir/$basename", "sample", "$sample");
            system(@cmd) && $self->throw("failed to run '@cmd'");
            $meta_count++;
        }

        $study = $vcf_meta->{study} if $vcf_meta->{study};
        if ($study) {
            @cmd = ("imeta", "add", "-d", "$idir/$basename", "study", "$study");
            system(@cmd) && $self->throw("failed to run '@cmd'");
            $meta_count++;
        }

        $release = $vcf_meta->{release} if $vcf_meta->{release};
        if ($release) {
            @cmd = ("imeta", "add", "-d", "$idir/$basename", "release", "$release");
            system(@cmd) && $self->throw("failed to run '@cmd'");
            $meta_count++;
        }

        # metadata check
        my $found_meta_count=0;
        my $pipe = "imeta ls -d $idir/$basename |";
        my $fh;
        open($fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
        while (<$fh>) {
            $found_meta_count++ if /attribute: /;
        }
        close($fh);
        $self->throw("Could not get all metadata for $idir/$basename") unless $meta_count == $found_meta_count;

        return 1;
    }

    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
}

1;
