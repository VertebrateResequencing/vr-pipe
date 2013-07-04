
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

class VRPipe::Steps::quantisnp_detect_cnv with VRPipe::StepRole {
    method options_definition {
        #sh run_quantisnp2.sh ./$QSNPDir/v79 --outdir $runQSNPDir/results --input-files $inFile --sampleid $shortFileName --levels ./$QSNPDir/config/levels.dat --config ./$QSNPDir/config/params.dat --chr 1:22 > $outFile
        
        # if optional given, not necessary, if not present, must be given
        
        return {
            run_quantisnp_script => VRPipe::StepOption->create(description => 'full path to run_quantisnp2.sh', optional => 1, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/run_quantisnp2.sh'),
            v79_dir              => VRPipe::StepOption->create(description => 'full path to v79 dir',           optional => 1, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/QSNP_dirs/v79'),
            levels_file          => VRPipe::StepOption->create(description => 'full path to levels.dat',        optional => 1, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/QSNP_dirs/config/levels.dat'),
            params_file          => VRPipe::StepOption->create(description => 'full path to params.dat',        optional => 1, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/QSNP_dirs/config/params.dat'),
            
            #out_dir => VRPipe::StepOption->create(description => 'full path to out dir', optional => 1, default_value => '/lustre/scratch102/user/pc12/genotyping/packages/QuantiSNP/kill_out/'),
            
            #reformat_annotation    => VRPipe::StepOption->create(description => 'Genome Studio annotation file (PC to expand this!)'),
            #reformat_mapping       => VRPipe::StepOption->create(description => 'file containing mapping of Genome Studio file columns to sample id'),
            #reformat_sample_number => VRPipe::StepOption->create(description => 'restrict reformatting to this number of samples',        optional => 1, default_value => 7),
        };
    }
    
    method inputs_definition {
        return {
            stepTwo_file_input_reformatted_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Reformatted Genome Studio file containing genotyping data'
            )
        };
    }
    
    #perl ./$PennCNVDir/detect_cnv.pl -test -hmm ./$PennCNVDir/lib/$hmmName -pfb ./$PennCNVDir/lib/$pfbFileName --confidence --region 1-22 --minsnp 1 -log $logFile -out $rawCNVFile $inFile > $outFile
    # sub = subroutine
    method body_sub {
        return sub {
            my $self                 = shift;
            my $options              = $self->options;
            my $run_quantisnp_script = $options->{run_quantisnp_script};
            my $v79_dir              = $options->{v79_dir};
            my $levels_file          = $options->{levels_file};
            my $params_file          = $options->{params_file};
            my $short_file_name      = '271298_A01_hipscigt5466706_50000';
            
            #my $out_dir                      = $options->{out_dir};
            
            my $quantisnp_detect_cnv_options = "--chr 1:22";
            #                            run_quantisnp2.sh ./$QSNPDir/v79 --outdir $runQSNPDir/results --input-files $inFile --sampleid $shortFileName --levels ./$QSNPDir/config/levels.dat --config ./$QSNPDir/config/params.dat --chr 1:22 > $outFile
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $reformatted_file (@{ $self->inputs->{stepTwo_file_input_reformatted_file} }) {
                my $reformatted_path = $reformatted_file->path;
                
                my $basename = $short_file_name . '.cnv';
                #my $out_dir       = $reformatted_file->dir;
                
                #print STDERR "#####BASENAME: $basename";
                #print STDERR "#####OUT_DIR:  $out_dir";
                
                my $quantisnp_cnv_file = $self->output_file(output_key => 'stepTwo_file_output_quantisnp_file', basename => "$basename", type => 'txt');
                #my $out_path      = $quantisnp_cnv_file->path;
                my $out_dir = $quantisnp_cnv_file->dir;
                
                #my $logFile       = $out_path . '.LOG';
                #my $cmd_line      = "perl $detect_cnv_script $detect_cnv_options -log $logFile -out $out_path $gs_path";
                
                #print STDERR "reformatted_file   = $reformatted_file\n";
                #print STDERR "reformatted_path   = $reformatted_path\n";
                #print STDERR "basename           = $basename\n";
                #print STDERR "quantisnp_cnv_file = $quantisnp_cnv_file\n";
                #print STDERR "out_path           = $out_path\n";
                #print STDERR "out_dir            = $out_dir\n";
                
                my $cmd_line = "sh $run_quantisnp_script $v79_dir --outdir $out_dir --input-files $reformatted_path --sampleid $short_file_name --levels $levels_file --config $params_file $quantisnp_detect_cnv_options";
                print STDERR "############ QUANTISNP DETECT_CNV: $cmd_line  #############################\n";
                $self->dispatch([$cmd_line, $req]); # RUN THE COMMAND
            }
        };
    }
    
    method outputs_definition {
        return {
            stepTwo_file_output_quantisnp_cnv_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'QuantiSNP CNV data',
                max_files   => -1,                  # -1 = As many as you like
                min_files   => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Detects CNVs using the Genome Studio genotyping files";
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
