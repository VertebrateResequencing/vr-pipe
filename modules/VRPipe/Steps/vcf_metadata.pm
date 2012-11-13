
=head1 NAME

VRPipe::Steps::vcf_metadata - a step

=head1 DESCRIPTION

Associates metadata with vcf files in the VRPipe database, making the file usable in other vcf-related Steps.
A basic implementation which may be mofified in future to include vcf-check output.

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

class VRPipe::Steps::vcf_metadata with VRPipe::StepRole {
    
    method options_definition {
        return {
            'md5sum_exe' => VRPipe::StepOption->create( description  => 'path to md5sum executable', optional => 1, default_value => 'md5sum'),
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
            my $md5sum_exe = $options->{md5sum_exe};
            my $cat_exe;
            
            my $req = $self->new_requirements(memory => 100, time => 1);

            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {

                my $input_path = $vcf_file->path;
                my $basename = $vcf_file->basename;
                my $cat_exe = $basename =~ /\.vcf.gz$/ ? 'zcat' : 'cat';

                my $cmd = "use VRPipe::Steps::vcf_metadata; VRPipe::Steps::vcf_metadata->add_metadata('$input_path', '$md5sum_exe', '$cat_exe');";
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
        return "Takes a vcf file and associates metadata with the file in the VRPipe database, making the vcf file usable in other vcf-related Steps";
    }
    
    method add_metadata (ClassName|Object $self: Str $input_vcf_path, Str $md5sum_exe, Str $cat_exe) {

        my $vcf_file = VRPipe::File->get(path => $input_vcf_path);

        # Add md5 checksum meta
        my $pipe = "$md5sum_exe $input_vcf_path |";
        open(my $fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
        my $md5_chksum;
        while (<$fh>) { # should only be one line
            ($md5_chksum,undef) = split;
        }
        close($fh);
        
        $self->throw("Could not get md5 checksum for $input_vcf_path") unless $md5_chksum;
        $self->throw("Got invalid md5 checksum $md5_chksum for $input_vcf_path") unless length($md5_chksum) == 32;

        $vcf_file->add_metadata({ md5 => $md5_chksum });
    
        # Get sample names from vcf header, store metadata as string of sample names delimited by '#'

        $pipe = "$cat_exe $input_vcf_path | head -500| grep '^#CHROM' |";
        open($fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
        my $samples;
        while (<$fh>) { # should only be one line
           /FORMAT\t(.*)$/;
            $samples = $1;
            $samples =~ s/\t/#/g;
        }
        close($fh);
        
        $self->throw("Could not get samples from $input_vcf_path") unless $samples;

        $vcf_file->add_metadata({ samples => $samples });

        return 1;
    }

    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
}

1;
