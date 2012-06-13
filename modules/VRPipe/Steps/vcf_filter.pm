=head1 NAME

VRPipe::Steps::vcf_filter - a step

=head1 DESCRIPTION

*** more documentation to come

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

class VRPipe::Steps::vcf_filter with VRPipe::StepRole {
	method options_definition {
		return { 'vcf-filter_file' => VRPipe::StepOption->get(description => 'vcf-filter parameters file'),
			'vcf-filter_exe' => VRPipe::StepOption->get(description => 'path to vcf-filter executable',
					optional => 1,
					default_value => 'vcf-filter'),
			'tabix_exe' => VRPipe::StepOption->get(description => 'path to your tabix executable',
					optional => 1,
					default_value => 'tabix') };
	}
	method inputs_definition {
		return { vcf_files => VRPipe::StepIODefinition->get(type => 'vcf',
				description => 'vcf files',
				max_files => -1) };
	}
	method body_sub {
		return sub {
			my $self = shift;

			my $options = $self->options;
			my $tabix_exe = $options->{tabix_exe};
			my $filter_exe = $options->{'vcf-filter_exe'};
			my $filter_file = $options->{'vcf-filter_file'};
			my $cat_exe;

			my $req = $self->new_requirements(memory => 500, time => 1);
			foreach my $vcf_file (@{$self->inputs->{vcf_files}}) {
				my $basename = $vcf_file->basename;
				if ($basename =~ /\.vcf.gz$/) {
					$basename =~ s/\.vcf.gz$/.filt.vcf.gz/;
					$cat_exe = 'zcat';
				}
				else {
					$basename =~ s/\.vcf$/.filt.vcf/;
					$cat_exe = 'cat';
				}

				my $filtered_vcf = $self->output_file(output_key => 'filtered_vcf', basename => $basename, type => 'vcf');
				my $tbi = $self->output_file(output_key => 'tbi_file', basename => $basename.'.tbi', type => 'bin');

				my $input_path = $vcf_file->path;
				my $output_path = $filtered_vcf->path;

				my $this_cmd = "$cat_exe $input_path | $filter_exe -f $filter_file | bgzip -c > $output_path; $tabix_exe -f -p vcf $output_path";

				$self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_filter', 'filter_vcf', [$this_cmd, $req, {output_files => [$filtered_vcf, $tbi]}]);
			}
		};
	}
	method outputs_definition {
		return { filtered_vcf => VRPipe::StepIODefinition->get(type => 'vcf',
				description => 'a filtered vcf file',
				max_files => -1),
			   tbi_file => VRPipe::StepIODefinition->get(type => 'bin',
					   description => 'a tbi file',
					   max_files => -1) };
	}
	method post_process_sub {
		return sub { return 1; };
	}
	method description {
		return "Soft-filters input VCFs, creating an output VCF for every input";
	}
	method max_simultaneous {
		return 0; # meaning unlimited
	}

	method filter_vcf (ClassName|Object $self: Str $cmd_line) {
		my ($input_path, $output_path) = $cmd_line =~ /^\S+ (\S+) .* vcf (\S+)$/;
		my $input_file = VRPipe::File->get(path => $input_path);

        my $input_recs = $input_file->num_records;

		$input_file->disconnect;
		system($cmd_line) && $self->throw("failed to run [$cmd_line]");

		my $output_file = VRPipe::File->get(path => $output_path);
		$output_file->update_stats_from_disc;
		my $output_recs = $output_file->num_records;

		unless ($output_recs == $input_recs) {
			$output_file->unlink;
			$self->throw("Output VCF has different number of data lines from input (input $input_recs, output $output_recs)");
		}
		else {
			return 1;
		}
	}
}

1;
