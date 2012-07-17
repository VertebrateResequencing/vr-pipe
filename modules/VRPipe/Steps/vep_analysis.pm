=head1 NAME

VRPipe::Steps::vep_analysis - a step

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

class VRPipe::Steps::vep_analysis with VRPipe::StepRole {
    method options_definition {
        return { 'vep_options' => VRPipe::StepOption->create(description => 'options to vep, excluding -i or -o'),
                 'vep_exe' => VRPipe::StepOption->create(description => 'path to your vep executable',
                                                               optional => 1,
                                                               default_value => 'variant_effect_predictor.pl') };
    }
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf',
                                                            description => 'vcf files',
                                                            max_files => -1) };
    }
	method body_sub {
		return sub {
			my $self = shift;

			my $options = $self->options;
			my $vep_exe = $options->{'vep_exe'};
			my $vep_opts = $options->{'vep_options'};
			my $cat_exe;

			if ($vep_opts =~ /-[i,o] /) {
				$self->throw("vep_options should not include the -i or -o option");
			}

			my $req = $self->new_requirements(memory => 5000, time => 1);
			foreach my $vcf_file (@{$self->inputs->{vcf_files}}) {
				my $basename = $vcf_file->basename;
				if ($basename =~ /\.vcf.gz$/) {
					$basename =~ s/\.vcf.gz$/.vep.txt/;
					$cat_exe = 'zcat';
				}
				else {
					$basename .= '.vep.txt';
					$cat_exe = 'cat';
				}
				my $vep_txt = $self->output_file(output_key => 'vep_txt', basename => $basename, type => 'txt',
													metadata => {source_vcf => $vcf_file->path->stringify});

				my $input_path = $vcf_file->path;
				my $output_path = $vep_txt->path;

				my $this_cmd = "$cat_exe $input_path | $vep_exe $vep_opts -o $output_path";

				$self->dispatch_wrapped_cmd('VRPipe::Steps::vep_analysis', 'vep_analysis', [$this_cmd, $req, {output_files => [$vep_txt]}]);
			}
		};
	}
    method outputs_definition {
        return { vep_txt => VRPipe::StepIODefinition->create(type => 'txt',
                                                             description => 'vep analysis output file',
                                                             max_files => -1) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Analyse VCF files with VEP";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method vep_analysis (ClassName|Object $self: Str $cmd_line) {

        my ($input_path, $output_path) = $cmd_line =~ /^\S+ (\S+) .* (\S+)$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        my $input_recs = $input_file->num_records;
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
		# Should at least be more than one output line per vcf record
        unless ($input_recs == 0 || $output_lines > $input_recs) {
            $output_file->unlink;
            $self->throw("VEP output has $output_lines lines, less than input vcf records $input_recs");
        }
        else {
            return 1;
        }
    }
}

1;
