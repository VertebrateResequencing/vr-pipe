
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

class VRPipe::Steps::quantisnp_reformat_gs_export with VRPipe::StepRole {
    method options_definition {
        # if optional given, not necessary, if not present, must be given
        
        return {
            reformat_script => VRPipe::StepOption->create(description => 'full path to convert_fcr_to_quanti-per_sample.pl',       optional => 0, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/convert_fcr_to_quanti-per_sample.pl'),
            manifest_file   => VRPipe::StepOption->create(description => 'full path to manifest file HumanCoreExome-12v1-0_A.csv', optional => 0, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/HumanCoreExome-12v1-0_A.csv'),
            
            #reformat_annotation    => VRPipe::StepOption->create(description => 'Genome Studio annotation file (PC to expand this!)'),
            #reformat_mapping       => VRPipe::StepOption->create(description => 'file containing mapping of Genome Studio file columns to sample id'),
            #reformat_sample_number => VRPipe::StepOption->create(description => 'restrict reformatting to this number of samples',        optional => 1, default_value => 7),
        };
    }
    
    method inputs_definition {
        return {
            stepOne_file_input_GS_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Genome Studio file containing genotyping data for reformatting'
            )
        };
    }
    
    #perl ./$PennCNVDir/detect_cnv.pl -test -hmm ./$PennCNVDir/lib/$hmmName -pfb ./$PennCNVDir/lib/$pfbFileName --confidence --region 1-22 --minsnp 1 -log $logFile -out $rawCNVFile $inFile > $outFile
    # sub = subroutine
    method body_sub {
        return sub {
            my $self            = shift;
            my $options         = $self->options;
            my $reformat_script = $options->{reformat_script};
            my $manifest_file   = $options->{manifest_file};
            #my $detect_cnv_pfb         = $options->{detect_cnv_pfb};
            #my $detect_cnv_options       = "-test -hmm $detect_cnv_hmm -pfb $detect_cnv_pfb --confidence --region 1-22 --minsnp 1";
            my $req = $self->new_requirements(memory => 500, time => 3);
            foreach my $gs_file (@{ $self->inputs->{stepOne_file_input_GS_file} }) {
                my $gs_path          = $gs_file->path;
                my $basename         = $gs_file->basename . '.reformat';
                my $reformatted_file = $self->output_file(output_key => 'stepOne_file_output_reformatted_file', basename => "$basename", type => 'txt');
                my $out_path         = $reformatted_file->path;
                #my $logFile       = $out_path . '.LOG';
                my $cmd_line = "perl $reformat_script $gs_path $manifest_file > $out_path";
                print STDERR "############ REFORMAT: $cmd_line  #############################\n";
                $self->dispatch([$cmd_line, $req]); # RUN THE COMMAND
            }
        };
    }
    
    method outputs_definition {
        return {
            stepOne_file_output_reformatted_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Reformatted file',
                max_files   => -1,                  # -1 = As many as you like
                min_files   => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Detects raw CNVs using the Genome Studio genotyping files";
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
