=head1 NAME

VRPipe::Steps::bam_sample_meta - a step

=head1 DESCRIPTION

Adds Sample and Gender meta-data to bam files

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::bam_sample_meta with VRPipe::StepRole {
    method options_definition {
        return { sample_sex_file => VRPipe::StepOption->get(description => 'File listing the sex (eg M or F) of samples, assumed_sex is used if no file provided', optional => 1),
                 assumed_sex => VRPipe::StepOption->get(description => 'If sex is not present for a sample in the sample sex file (or no file provided), then this sex is assumed', optional => 1, default_value => 'F'),
        };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->get(type => 'bam', max_files => -1, description => '1 or more bam files') };
    }
	method body_sub {
		return sub {
			my $self = shift;

			my $options = $self->options;
			my $assumed_sex = $options->{assumed_sex};
			my $sample_sex_file = $options->{sample_sex_file};

			my $sample_sex_path='';
			if ($sample_sex_file) {
				$sample_sex_path = Path::Class::File->new($sample_sex_file);
				$self->throw("sample_sex_file must be an absolute path") unless $sample_sex_path->is_absolute;
			}

            my $req = $self->new_requirements(memory => 500, time => 1);

			foreach my $bamfile (@{$self->inputs->{bam_files}}) {
				my $bam_path = $bamfile->path;
				my $cmd = "use VRPipe::Steps::bam_sample_meta; VRPipe::Steps::bam_sample_meta->add_sample_meta('$bam_path', '$sample_sex_path', '$assumed_sex');";
				$self->dispatch_vrpipecode($cmd, $req);

			}
		};
	}
    method outputs_definition {
        return {
        };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Adds Sample and Gender meta-data to bam files";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }

	method add_sample_meta (ClassName|Object $self: Str|File $bam!, Str|File $sample_sex_path, Str $assumed_sex) {

		my $bam_sample;
		open(my $bamhdrs,"samtools view -H $bam |") or $self->throw("samtools view -H $bam: $!");
		while (my $bamhdr = <$bamhdrs>) {
			next unless $bamhdr=~/^\@RG/;
			if ( $bamhdr=~/SM:(\S+)/ ) {
				$bam_sample = $1;
				last; # Assuming single-sample bam
			}
		}
		close($bamhdrs);

		my $bam_sex;
		if ($sample_sex_path) {
			my $sex_file = VRPipe::File->get(path => $sample_sex_path);
			my $fh = $sex_file->openr;
			while (<$fh>) {
				chomp;
				my ($sample, $sex) = split /\t/;
				next unless $sample eq $bam_sample;
				$bam_sex = $sex;
			}
			$fh->close;
		}
		$bam_sex ||= $assumed_sex;

		my $bamfile = VRPipe::File->get(path => $bam);
		$bamfile->add_metadata({sample=> $bam_sample, sex => $bam_sex});

	}
}

1;
