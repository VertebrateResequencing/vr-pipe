
=head1 NAME

VRPipe::Steps::split_genome_studio_genotype_files - a step

=head1 DESCRIPTION

This step obtains the metadata from gtc files that have previously been
imported from iRODS. The sample metadata is utilised in a zgrep command to get
sample-specific genotype data from a gzipped file of genome studio genotyping
data of multiple samples involved in the current analysis.

The output file contains the genotype data for the particular sample and has
the sample metadata attached to it.

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::split_genome_studio_genotype_files  with VRPipe::StepRole  {
    method options_definition {
        return {
            zgrep_exe        => VRPipe::StepOption->create(description => 'full path to zgrep command that is used to retrieve the header and sample data from the gzipped genome studio genotyping data file', optional => 1, default_value => 'zgrep'),
            header_regex     => VRPipe::StepOption->create(description => 'regex of a field in the header line of the gzipped genome studio genotype file that is used to retrieve the header using zgrep',     optional => 1, default_value => 'Allele1'),
            reheader_penncnv => VRPipe::StepOption->create(description => 'optionally, header the genotype file with a penncnv-friendly format',                                                                optional => 1)
        };
    }
    
    method inputs_definition {
        return {
            gtc_files => VRPipe::StepIODefinition->create(
                type            => 'gtc',
                description     => 'gtc file with associated genome studio irods metadata, where analysis files have already been downloaded',
                max_files       => -1,
                check_existence => 0,
                metadata        => { sample => 'sample name for cell line', infinium_sample => 'the sample id found in the corresponding fcr file noted in irods_analysis_files', irods_analysis_files => 'full irods path to multi-sample genome studio analysis fcr file', irods_local_storage_dir => 'local base directory where the irods_analysis_files were downloaded' }
            ),
        };
    }
    
    method outputs_definition {
        return {
            fcr_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a single-sample fcr file',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self             = shift;
            my $options          = $self->options;
            my $zgrep_exe        = $options->{zgrep_exe};
            my $header_regex     = $options->{header_regex};
            my $reheader_penncnv = $options->{reheader_penncnv} ? $options->{reheader_penncnv} : undef;
            my $req              = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $gtc_file (@{ $self->inputs->{gtc_files} }) {
                my $meta           = $gtc_file->metadata;
                my $sample         = $meta->{sample};
                my $fcr_sample     = $meta->{infinium_sample};
                my $analysis_files = $meta->{irods_analysis_files};
                my $fcr_file;
                if (ref($analysis_files)) {
                    foreach my $file (@$analysis_files) {
                        if ($file =~ /\.fcr/i) {
                            $fcr_file = $file;
                            last;
                        }
                    }
                }
                else {
                    $fcr_file = $analysis_files;
                }
                my $local_dir             = $meta->{irods_local_storage_dir};
                my $multi_sample_fcr_file = file($local_dir, $fcr_file)->stringify;
                my $this_grep_exe         = $multi_sample_fcr_file =~ /\.gz$/ ? $zgrep_exe : 'grep';
                my $basename              = $sample . '.genotyping.fcr.txt';
                my $sample_genotype_file  = $self->output_file(output_key => 'fcr_files', basename => $basename, type => 'txt', metadata => $meta)->path;
                my $header_cmd            = $reheader_penncnv ? "$zgrep_exe $header_regex $reheader_penncnv > $sample_genotype_file " : "$this_grep_exe $header_regex $multi_sample_fcr_file > $sample_genotype_file ";
                my $cmd_line              = $header_cmd . qq[&& $this_grep_exe "$fcr_sample\\s" $multi_sample_fcr_file >> $sample_genotype_file];
                $self->dispatch([$cmd_line, $req]);
            }
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method description {
        return "This step uses zgrep to create sample genotype data text files from a gzipped genotype data file produced by genome studio that contains genotyping data from multiple samples.";
    }

}

1;
