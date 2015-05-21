
=head1 NAME

VRPipe::Steps::pluritest_reformat_genome_studio_expression_files - a step

=head1 DESCRIPTION

Converts the Genome Studio csv files into a format that is suitable for
processing by the PluriTest R package to determine pluripotency

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>. Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013,2014 Genome Research Limited.

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

class VRPipe::Steps::pluritest_reformat_genome_studio_expression_files with VRPipe::StepRole {
    method options_definition {
        return {
            reformat_exe           => VRPipe::StepOption->create(description => 'full path to genome_studio_gene_expression_reformat.pl', optional => 1, default_value => 'genome_studio_gene_expression_profile_reformat.pl'),
            reformat_sample_number => VRPipe::StepOption->create(description => 'restrict reformatting to this number of samples',        optional => 1, default_value => 7),
        };
    }
    
    method inputs_definition {
        return {
            annotation_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'the annotation file that is used by pluritest',
            ),
            mapping_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that maps the samples to the columns in the sample profile file',
            ),
            profile_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that contains the GenomeStudio profile for the samples',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                   = shift;
            my $options                = $self->options;
            my $reformat_exe           = $options->{reformat_exe};
            my $reformat_sample_number = $options->{reformat_sample_number};
            
            my ($profile_file)    = @{ $self->inputs->{profile_files} };
            my $profile_path      = $profile_file->path;
            my ($annotation_file) = @{ $self->inputs->{annotation_files} };
            my ($mapping_file)    = @{ $self->inputs->{mapping_files} };
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $reformat_options = "--annot " . $annotation_file->path . " --mapping " . $mapping_file->path . " --samples $reformat_sample_number";
            
            my $basename      = $profile_file->basename . '.reformat';
            my $reformat_file = $self->output_file(output_key => 'reformat_files', basename => "$basename", type => 'txt', metadata => $profile_file->metadata);
            my $out_path      = $reformat_file->path;
            
            my $cmd_line = "$reformat_exe --profile $profile_path $reformat_options --out $out_path";
            $self->dispatch_wrapped_cmd('VRPipe::Steps::pluritest_reformat_genome_studio_expression_files', 'reformat_gs_file', [$cmd_line, $req]);
        
        };
    }
    
    method outputs_definition {
        return {
            reformat_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Files with reformatted Genome Studio gene expression data',
                max_files   => -1,
                min_files   => 0,
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
        
        #reformatted file should have header line with ProbeID as the first item
        my $reformf    = $output_file->openr;
        my $first_line = <$reformf>;
        unless ($first_line =~ /^ProbeID/) {
            $output_file->unlink;
            $self->throw("Reformatted file does not have correct header with Probe ID as the first expected column\n");
        }
        # Should be one output line per Genome Studio record input
        if ($output_lines == $input_recs) {
            $self->relate_input_to_output($input_path, 'reformat_for_pluritest', $output_path);
            return 1;
        }
        else {
            $output_file->unlink;
            $self->throw("The reformatted file does not have the same number of records as the Genome Studio input file");
        }
    }
}

1;
