=head1 NAME

VRPipe::Steps::sex_to_ploidy - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::sex_to_ploidy with VRPipe::StepRole {
    method options_definition {
        return { sample_sex_file => VRPipe::StepOption->create(description => 'File listing the sex (M or F) of samples'),
                 assumed_sex => VRPipe::StepOption->create(description => 'If M or F is not present for a sample in the sample sex file, then this sex is assumed', optional => 1, default_value => 'F') };
    }
    method inputs_definition {
        return { chunked_regions_file => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'Chromosomal regions file split according to chunk size'),
                 bcf_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'The bcf files which were chunked by chunked_regions file'), };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $assumed_sex = $options->{assumed_sex};
            my $sample_sex_file = Path::Class::File->new($options->{sample_sex_file});
            $self->throw("sample_sex_file must be an absolute path") unless $sample_sex_file->is_absolute;
            
            my %bcfs;
            foreach my $bcf (@{$self->inputs->{bcf_files}}) {
                $bcfs{$bcf->metadata->{region}} = $bcf;
            }
            
            my $chunk_file = $self->inputs->{chunked_regions_file}[0];
            my $fh = $chunk_file->openr;
            while (<$fh>) {
                my ($chr, $from, $to, $female_ploidy, $male_ploidy) = split /\t/;
                my $region = "${chr}_${from}-${to}";
                next unless (exists $bcfs{$region});
                my $bcf = $bcfs{$region};
                my $base = $bcf->basename;
                $base =~ s/bcf$/samples/;
                my $sample_ploidy_file = $self->output_file(output_key => 'sample_ploidy_files',
                                                            basename => $base,
                                                            type => 'txt', 
                                                            metadata => { source_bcf => $bcf->path->stringify });
                my $sample_ploidy_path = $sample_ploidy_file->path;
                my $req = $self->new_requirements(memory => 50, time => 1);
                my $cmd = "use VRPipe::Steps::sex_to_ploidy; VRPipe::Steps::sex_to_ploidy->write_sample_ploidy_file('$sample_sex_file', '$sample_ploidy_path', female_ploidy => $female_ploidy, male_ploidy => $male_ploidy, assumed_sex => '$assumed_sex');";
                $self->dispatch_vrpipecode($cmd, $req, {output_files => [$sample_ploidy_file]});
            }
        };
    }
    method outputs_definition {
        return { sample_ploidy_files => VRPipe::StepIODefinition->create(type => 'txt', 
                                                                      description => 'Files listing sample name and ploidy for each region to be used as input to bcftools calling', 
                                                                      max_files => -1,
                                                                      metadata => { source_bcf => 'bcf for which this samples file was generated' }) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Creates files listing sample and ploidy to be used as input to bcftools calling";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }

    method write_sample_ploidy_file (ClassName|Object $self: Str|File $sample_sex_list, Str|File $sample_ploidy_path, Str :$female_ploidy, Str :$male_ploidy, Str :$assumed_sex) {
        
        my $sex_file = VRPipe::File->create(path => $sample_sex_list);
        my $fh = $sex_file->openr;
        
        my $ploidy_file = VRPipe::File->get(path => $sample_ploidy_path);
        my $pfh = $ploidy_file->openw;
        
        $ploidy_file->disconnect;
        my $expected_samples = 0;
        while (<$fh>) {
            chomp;
            my ($sample, $sex) = split /\t/;
            $sex ||= $assumed_sex;
            my $ploidy = $sex eq 'F' ? $female_ploidy : $male_ploidy;
            next unless $ploidy;
            print $pfh "$sample\t$ploidy\n";
            $expected_samples++;
        }
        $pfh->close;
        $ploidy_file->update_stats_from_disc(retries => 3);
        
        # check number of samples is as expected
        my $actual_samples = $ploidy_file->lines;
        if ($actual_samples == $expected_samples) {
            return 1;
        } else {
            # $ploidy_file->unlink;
            $self->throw("write_sample_ploidy_file failed because $actual_samples samples were generated in the output ploidy files, yet we expected $expected_samples samples");
        }
    }
}

1;
