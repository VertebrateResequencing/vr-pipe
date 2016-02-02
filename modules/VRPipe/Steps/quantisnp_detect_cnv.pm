
=head1 NAME

VRPipe::Steps::quantisnp_detect_cnv - a step

=head1 DESCRIPTION

Detects CNVs using reformatted Genome Studio genotyping fcr files

=head1 AUTHOR

Phil Carter <pc12@sanger.ac.uk>, John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Steps::quantisnp_detect_cnv with VRPipe::StepRole {
    method options_definition {
        return {
            run_quantisnp_script => VRPipe::StepOption->create(description => 'full path to run_quantisnp2.sh', optional => 0),
            quantisnp_options    => VRPipe::StepOption->create(description => 'options to run_quantisnp2.sh',   optional => 0, default => "--chr 1:22"),
            v79_dir              => VRPipe::StepOption->create(description => 'full path to v79 dir',           optional => 0),
            levels_file          => VRPipe::StepOption->create(description => 'full path to levels.dat',        optional => 0),
            params_file          => VRPipe::StepOption->create(description => 'full path to params.dat',        optional => 0),
        };
    }
    
    method inputs_definition {
        return {
            fcr_files_reformatted => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Reformatted Genome Studio fcr files containing genotyping data'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self                         = shift;
            my $options                      = $self->options;
            my $run_quantisnp_script         = $options->{run_quantisnp_script};
            my $v79_dir                      = $options->{v79_dir};
            my $levels_file                  = $options->{levels_file};
            my $params_file                  = $options->{params_file};
            my $quantisnp_detect_cnv_options = $options->{quantisnp_options};
            my $req                          = $self->new_requirements(memory => 500, time => 1);
            
            foreach my $reformatted_file (@{ $self->inputs->{fcr_files_reformatted} }) {
                my $reformatted_path   = $reformatted_file->path;
                my $reformat_meta      = $reformatted_file->metadata;
                my $sample             = $reformat_meta->{sample};
                my $basename           = $sample . '.cnv';
                my $quantisnp_cnv_file = $self->output_file(output_key => 'stepTwo_file_output_quantisnp_file', basename => "$basename", type => 'txt', metadata => $reformat_meta);
                my $out_dir            = $quantisnp_cnv_file->dir;
                my $cmd_line           = "sh $run_quantisnp_script $v79_dir --outdir $out_dir --input-files $reformatted_path --sampleid $sample --levels $levels_file --config $params_file $quantisnp_detect_cnv_options";
                $self->dispatch([$cmd_line, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            quantisnp_cnv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'QuantiSNP CNV data',
                max_files   => -1,
                min_files   => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Detects CNVs using reformatted Genome Studio fcr files";
    }
    
    method max_simultaneous {
        return 0;
    }

}

1;
