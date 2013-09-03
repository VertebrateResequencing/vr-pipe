
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
            zgrep_exe          => VRPipe::StepOption->create(description => 'full path to zgrep command that is used to retrieve the header and sample data from the gzipped genome studio genotyping data file', optional => 1, default_value => 'zgrep'),
            header_regex       => VRPipe::StepOption->create(description => 'regex of a field in the header line of the gzipped genome studio genotype file that is used to retrieve the header using zgrep',     optional => 1, default_value => 'Allele1'),
            reheader_penncnv   => VRPipe::StepOption->create(description => 'optionally, header the genotype file with a penncnv-friendly format',                                                                optional => 1),
            external_gzip_file => VRPipe::StepOption->create(description => 'optionally, provide path to a gzipped genotype file that can override the path provided in the input gtc file metadata',             optional => 1),
        };
    }
    
    method inputs_definition {
        return {
            gtc_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'gtc file with associated genome studio irods file path and analysis_uuid metadata',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' }
            ),
        };
    }
    
    method outputs_definition {
        return {
            gtype_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file on a local disc, extracted from iRODs',
                max_files   => -1,
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self               = shift;
            my $options            = $self->options;
            my $zgrep_exe          = $options->{zgrep_exe};
            my $header_regex       = $options->{header_regex};
            my $reheader_penncnv   = $options->{reheader_penncnv} ? $options->{reheader_penncnv} : undef;
            my $external_gzip_file = $options->{external_gzip_file} ? $options->{external_gzip_file} : undef;
            my $req                = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $gtc_file (@{ $self->inputs->{gtc_files} }) {
                my $sample     = $gtc_file->metadata->{sample};
                my $library    = $gtc_file->metadata->{library};
                my $lib_sample = (split /_/, $library)[-1];
                #the gzipped file has already been downloaded by the iRODS vrtrack updater and the path is stored in the storage_path metadata field
                #alternatively, it can be overridden by an external gzipped genotype file if required
                my $genotype_gzip_path   = $external_gzip_file ? $external_gzip_file : $gtc_file->metadata->{storage_path};
                my $basename             = $sample . '.genotyping.fcr.txt';
                my $sample_genotype_file = $self->output_file(output_key => 'gtype_files', basename => $basename, type => 'txt', metadata => $gtc_file->metadata)->path;
                my $header_cmd           = $reheader_penncnv ? "$zgrep_exe $header_regex $reheader_penncnv > $sample_genotype_file " : "$zgrep_exe $header_regex $genotype_gzip_path > $sample_genotype_file ";
                my $cmd_line             = $header_cmd . "&& $zgrep_exe $lib_sample $genotype_gzip_path >> $sample_genotype_file";
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
