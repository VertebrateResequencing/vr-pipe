
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

class VRPipe::Steps::genome_studio_expression_reformat with VRPipe::StepRole {
    method options_definition {
        return {
            reformat_exe           => VRPipe::StepOption->create(description => 'full path to genome_studio_gene_expression_reformat.pl', optional => 1, default_value => 'genome_studio_gene_expression_reformat.pl'),
            reformat_annotation    => VRPipe::StepOption->create(description => 'Genome Studio annotation file (PC to expand this!)'),
            reformat_mapping       => VRPipe::StepOption->create(description => 'file containing mapping of Genome Studio file columns to sample id'),
            reformat_sample_number => VRPipe::StepOption->create(description => 'restrict reformatting to this number of samples',        optional => 1, default_value => 7),
        };
    }
    
    method inputs_definition {
        return {
            gs_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Genome Studio file containing gene expression data that needs to be reformatted'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self                   = shift;
            my $options                = $self->options;
            my $reformat_exe           = $options->{reformat_exe};
            my $reformat_annotation    = $options->{reformat_annotation};
            my $reformat_mapping       = $options->{reformat_mapping};
            my $reformat_sample_number = $options->{reformat_sample_number};
            my $reformat_options       = "--annot $reformat_annotation --mapping $reformat_mapping --samples $reformat_sample_number";
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $gs_file (@{ $self->inputs->{gs_file} }) {
                my $gs_path       = $gs_file->path;
                my $basename      = $gs_file->basename . '.reformat';
                my $reformat_file = $self->output_file(output_key => 'reformat_files', basename => "$basename", type => 'txt');
                my $out_path      = $reformat_file->path;
                my @output_files  = ($reformat_file);
                my $cmd_line      = "$reformat_exe --profile $gs_path $reformat_options --out $out_path";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::genome_studio_expression_reformat', 'reformat_gs_file', [$cmd_line, $req, { output_files => [$reformat_file] }]);
            }
        
        };
    }
    
    method outputs_definition {
        return {
            reformat_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Files with reformatted Genome Studio gene expression data',
                max_files   => -1,
                min_files   => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts the Genome Studio gene expression files into a format that is suitable for processing by the PluriTest R package";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method reformat_gs_file (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /--profile (\S+) .* --out (\S+)$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        my $input_recs = $input_file->num_records;
        $input_file->disconnect;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        my $reformf    = $output_file->openr;
        my $first_line = <$reformf>;
        #reformatted file should have header with ProbeID as the first item
        unless ($first_line =~ /^ProbeID/) {
            $output_file->unlink;
            $self->throw("Reformatted file does not have correct header with Probe ID as the first expected column\n");
        }
        # Should be one output line per Genome Studio record input
        if ($output_lines == $input_recs) {
            return 1;
        }
        else {
            $output_file->unlink;
            $self->throw("The reformatted file does not have the same number of records as the Genome Studio input file");
        }
    }
}

1;
