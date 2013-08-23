
=head1 NAME

VRPipe::Steps::pluritest_reformat_genome_studio_expression_files - a step

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
                description => 'the annotation file that is used by pluritest - if more than one is provided, the intersection of the annotation files is produced ',
                max_files   => -1
            ),
            mapping_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that maps the samples to the columns in the sample profile file',
                max_files   => -1
            ),
            profile_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a file that contains the GenomeStudio profile for the samples',
                max_files   => -1,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self                   = shift;
            my $options                = $self->options;
            my $reformat_exe           = $options->{reformat_exe};
            my $reformat_sample_number = $options->{reformat_sample_number};
            
            # group the profile with its mapping and annotation files
            my %by_tag;
            foreach my $file (@{ $self->inputs->{profile_files} }, @{ $self->inputs->{annotation_files} }, @{ $self->inputs->{mapping_files} }) {
                push(@{ $by_tag{ $file->metadata->{merge_tag_id} } }, $file->path);
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            while (my ($tag, $files) = each %by_tag) {
                $self->throw("There was not exactly 1 annotation file and 1 mapping file and 1 profile file per merge tag id $tag (@$files)") unless @$files == 3;
                my $reformat_options = "--annot " . $files->[1] . " --mapping " . $files->[2] . " --samples $reformat_sample_number";
                my $profile_path     = $files->[0];
                my $profile_file     = VRPipe::File->get(path => $profile_path);
                my $basename         = $profile_file->basename . '.reformat';
                my $reformat_file    = $self->output_file(output_key => 'reformat_files', basename => "$basename", type => 'txt', metadata => $profile_file->metadata);
                my $out_path         = $reformat_file->path;
                my $cmd_line         = "$reformat_exe --profile $profile_path $reformat_options --out $out_path";
                $self->dispatch_wrapped_cmd('VRPipe::Steps::pluritest_reformat_genome_studio_expression_files', 'reformat_gs_file', [$cmd_line, $req, { output_files => [$reformat_file] }]);
            }
        
        };
    }
    
    method outputs_definition {
        return {
            reformat_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Files with reformatted Genome Studio gene expression data',
                max_files   => -1,
                min_files   => 0,
                metadata    => { merge_tag_id => 'tag id to enable sample to be identified in the multi-sample profile file', lanes => 'comma-separated list of lanes that the pluritest analysis is being performed on' },
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
            return 1;
        }
        else {
            $output_file->unlink;
            $self->throw("The reformatted file does not have the same number of records as the Genome Studio input file");
        }
    }
}

1;
