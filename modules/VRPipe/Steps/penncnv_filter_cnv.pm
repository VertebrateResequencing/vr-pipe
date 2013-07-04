
=head1 NAME

VRPipe::Steps::genome_studio_expression_reformat - a step

=head1 DESCRIPTION

Converts the Genome Studio csv files into a format that is suitable for
processing by the PluriTest R package to determine pluripotency

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

class VRPipe::Steps::penncnv_filter_cnv with VRPipe::StepRole {
    method options_definition {
        # if optional given, not necessary, if not present, must be given
        
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
    
    #filter_cnv.pl -numsnp 10 -length 130k rawcnv//271298_G06_HAPMAP5265757-Exome1.1-custom-42.rawcnv --confidence 10
    
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
                print STDERR "############ FILTER_CNV: $cmd_line  #############################\n";
                $self->dispatch([$cmd_line, $req]); # RUN THE COMMAND
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
    
    # actually running command
    
    #    method detect_cnv (ClassName|Object $self: Str $cmd_line) {
    #        #my ($input_path, $output_path) = $cmd_line =~ /--profile (\S+) .* --out (\S+)$/;
    #        #my $input_file = VRPipe::File->get(path => $input_path);
    #        #my $input_recs = $input_file->num_records;
    #        #$input_file->disconnect; # finished with the file, kind of a close but not really
    #
    #        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
    #
    #        my $output_file = VRPipe::File->get(path => $output_path);
    #        $output_file->update_stats_from_disc; # get some info so you can pull stuff like lines, # records, etc
    #        #my $output_lines = $output_file->lines;
    #
    #        #reformatted file should have header line with ProbeID as the first item
    #        #my $reformf    = $output_file->openr; # standard open within vrpipe
    #        #my $first_line = <$reformf>;
    #        #unless ($first_line =~ /^ProbeID/) {
    #           $output_file->unlink; # kind of disconnecting
    #        #    $self->throw("Reformatted file does not have correct header with Probe ID as the first expected column\n");
    #        #}
    #        # Should be one output line per Genome Studio record input
    #        #if ($output_lines == $input_recs) {
    #        #    return 1;
    #        #}
    #        #else {
    #        #    $output_file->unlink;
    #        #    $self->throw("The reformatted file does not have the same number of records as the Genome Studio input file");
    #        #}
    #    }
}

1;
