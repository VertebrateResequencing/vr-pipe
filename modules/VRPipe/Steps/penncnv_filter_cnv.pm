
=head1 NAME

VRPipe::Steps::penncnv_filter_cnv - a step

=head1 DESCRIPTION

Filters raw CNVs found by PennCNV detect_cnv.pl

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

class VRPipe::Steps::penncnv_filter_cnv with VRPipe::StepRole {
    method options_definition {
        return {
            filter_cnv_script => VRPipe::StepOption->create(description => 'full path to filter_cnv.pl', optional => 1, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/PennCNV/PennCNV/filter_cnv.pl'),
        };
    }
    
    method inputs_definition {
        return {
            stepTwo_file_input_raw_cnv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Raw CNV file produced by detect_cnv.pl within PennCNV',
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self               = shift;
            my $options            = $self->options;
            my $filter_cnv_script  = $options->{filter_cnv_script};
            my $filter_cnv_options = "--numsnp 10 --length 130k --confidence 10";
            my $req                = $self->new_requirements(memory => 500, time => 1);
            foreach my $raw_cnv_file (@{ $self->inputs->{stepTwo_file_input_raw_cnv_file} }) {
                my $raw_cnv_path    = $raw_cnv_file->path;
                my $basename        = $raw_cnv_file->basename . '.filtercnv';
                my $filter_cnv_file = $self->output_file(output_key => 'stepTwo_file_output_filter_cnv_file', basename => "$basename", type => 'txt', metadata => $raw_cnv_file->metadata);
                my $out_path        = $filter_cnv_file->path;
                my $cmd_line        = "perl $filter_cnv_script $filter_cnv_options --out $out_path $raw_cnv_path";
                $self->dispatch([$cmd_line, $req]);
            }
        };
    }
    
    method outputs_definition {
        return {
            stepTwo_file_output_filter_cnv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Filter CNV data',
                max_files   => -1,                                                                                                                    # -1 = As many as you like
                min_files   => 0,
                metadata    => { sample => 'sample name for cell line', storage_path => 'full path to iRODS file', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Filters raw CNVs found by PennCNV detect_cnv.pl";
    }
    
    method max_simultaneous {
        return 0;
    }

}

1;
