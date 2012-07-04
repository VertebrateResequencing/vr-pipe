=head1 NAME

VRPipe::Steps::vcf_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::vcf_index with VRPipe::StepRole {
    method options_definition {
        return {
		 'tabix_exe' => VRPipe::StepOption->create(description => 'path to tabix executable', optional => 1, default_value => 'tabix')
		};
    }
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', description => 'vcf files to index', max_files => -1) 
		};
    }
	method body_sub {
		return sub {
			my $self = shift;

			my $options = $self->options;
			my $tabix_exe = $options->{tabix_exe};

			my $req = $self->new_requirements(memory => 500, time => 1);
			foreach my $vcf (@{$self->inputs->{vcf_files}}) {

				my $vcf_path = $vcf->path;
				my $basename = $vcf->basename;
				my $tbi = $self->output_file(output_key => 'tbi_file', output_dir => $vcf->dir, basename => $basename.'.tbi', type => 'bin');

				my $cmd = "$tabix_exe -f -p vcf $vcf_path;";
				$self->dispatch([$cmd, $req, {output_files => [$tbi]}]); 
			}
		};
	}
    method outputs_definition {
        return {
                 tbi_file => VRPipe::StepIODefinition->create(type => 'bin', description => 'a tbi file', max_files => -1)
		};
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "tabix index VCF files";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
