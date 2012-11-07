
=head1 NAME

VRPipe::Steps::bin2hapmap_sites - a step

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

class VRPipe::Steps::bin2hapmap_sites with VRPipe::StepRole {
    method options_definition {
        return {
            hapmap2bin_sample_genotypes_file => VRPipe::StepOption->create(description => 'absolute path to binary file of sample genotypes, produced by hapmap2bin'),
            bin2hapmap_exe                   => VRPipe::StepOption->create(
                description   => 'path to bin2hapmap executable',
                optional      => 1,
                default_value => 'bin2hapmap'
            )
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self           = shift;
            my $options        = $self->options;
            my $gtype_snps_bin = file($options->{hapmap2bin_sample_genotypes_file});
            $self->throw("hapmap2bin_sample_genotypes_file must be an absolute path") unless $gtype_snps_bin->is_absolute;
            my $bin2hapmap  = $options->{bin2hapmap_exe};
            my $req         = $self->new_requirements(memory => 3900, time => 1);
            my $hapmap_file = $self->output_file(
                output_key => 'hapmap_file',
                output_dir => $gtype_snps_bin->dir->stringify,
                basename   => $gtype_snps_bin->basename . '.hapmap.txt',
                type       => 'txt'
            );
            my $hapmap_path = $hapmap_file->path;
            my $cmd         = qq[$bin2hapmap -l $gtype_snps_bin > $hapmap_path];
            $self->dispatch([$cmd, $req, { output_files => [$hapmap_file], block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return {
            hapmap_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'text file with the hapmap sites output from the snp binary file'
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Convert a hapmap2bin-generated snp binary file to a textual list of sites with bin2hapmap.";
    }
    
    method max_simultaneous {
        return 0;
    }
}

1;
