
=head1 NAME

VRPipe::Steps::cnv_control_comparison - a step

=head1 DESCRIPTION

calculates the intersection of the control bed file with those of the stem
cells.

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

class VRPipe::Steps::cnv_control_comparison with VRPipe::StepRole {
    method options_definition {
        return {
            bedtools_exe     => VRPipe::StepOption->create(description => 'path to your bedtools executable',                                                               optional => 1, default_value => 'bedtools'),
            cell_control_tag => VRPipe::StepOption->create(description => 'Tag assigned to the file metadata control parameter to determine that a cell line is a control', optional => 1, default_value => 'Control'),
        };
    }
    
    method inputs_definition {
        return {
            bed_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'Reformatted bed file from CNV analyses of stem cell samples (non-control)',
                metadata    => { sample => 'sample name for cell line', control => 'determine if cell line is control or stem cell', individual => 'cohort id', analysis_uuid => 'analysis_uuid' },
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self             = shift;
            my $options          = $self->options;
            my $bedtools_exe     = $options->{bedtools_exe};
            my $cell_control_tag = $options->{cell_control_tag};
            my $req              = $self->new_requirements(memory => 500, time => 1);
            my %control_files;
            my %stem_files;
            foreach my $bed_file (@{ $self->inputs->{bed_files} }) {
                my $bed_path   = $bed_file->path;
                my $meta       = $bed_file->metadata;
                my $individual = $meta->{individual};
                my $control    = $meta->{control};
                if ($control eq $cell_control_tag) {
                    $control_files{$individual} = $bed_path;
                }
                else {
                    push(@{ $stem_files{$individual} }, $bed_path);
                }
            }
            
            #check stem cells cohorts === control cohorts
            foreach my $cohort (keys %stem_files) {
                my $control_path    = $control_files{$cohort};
                my @stem_cell_paths = @{ $stem_files{$cohort} };
                foreach my $stem_path (@stem_cell_paths) {
                    my $stem_file          = VRPipe::File->get(path => $stem_path);
                    my $intersect_basename = $stem_file->basename . '.intersection';
                    my $diff_basename      = $stem_file->basename . '.diff';
                    my $intersect_file     = $self->output_file(output_key => 'intersect_file', basename => "$intersect_basename", type => 'txt', metadata => $stem_file->metadata);
                    my $intersect_path     = $intersect_file->path;
                    my $diff_file          = $self->output_file(output_key => 'diff_file', basename => "$diff_basename", type => 'txt', metadata => $stem_file->metadata);
                    my $diff_path          = $diff_file->path;
                    my $this_cmd           = "$bedtools_exe intersect -u -f 0.5 -r -a $stem_path -b $control_path > $intersect_path && $bedtools_exe intersect -v -f 0.5 -r -a $stem_path -b $control_path > $diff_path";
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::cnv_control_comparison', 'bed_intersection', [$this_cmd, $req, { output_files => [$intersect_file, $diff_file] }]);
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            intersect_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Bed file containing the intersection of control and stem cell bed files',
                max_files   => -1,                                                                       # -1 = As many as you like
                min_files   => 0
            ),
            diff_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'Bed file containing the intersection of control and stem cell bed files',
                max_files   => -1,                                                                       # -1 = As many as you like
                min_files   => 0
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Calculates the intersection between bed control and stem cell CNV ";
    }
    
    method max_simultaneous {
        return 0;
    }
    
    method bed_intersection (ClassName|Object $self: Str $cmd_line) {
        my ($cnv_stem_path, $intersection_path, $diff_path) = $cmd_line =~ /^.* -a (\S+) .* > (\S+) .* (\S+)$/;
        
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $cnv_stem_file = VRPipe::File->get(path => $cnv_stem_path);
        my $cnv_stem_lines = $cnv_stem_file->lines;
        $cnv_stem_file->disconnect;
        
        my $intersection_file = VRPipe::File->get(path => $intersection_path);
        $intersection_file->update_stats_from_disc;
        my $intersection_lines = $intersection_file->lines;
        
        my $diff_file = VRPipe::File->get(path => $diff_path);
        $diff_file->update_stats_from_disc;
        my $diff_lines = $diff_file->lines;
        
        my $combined_lines = $intersection_lines + $diff_lines;
        
        if ($combined_lines == $cnv_stem_lines) {
            my $new_meta = { cnv_minus_control => $diff_lines, cnv_total => $combined_lines };
            $intersection_file->add_metadata($new_meta);
            $diff_file->add_metadata($new_meta);
            $cnv_stem_file->add_metadata($new_meta);
            return 1;
        }
        else {
            $self->warn("Number of intersection and diff file lines do not equal the original stem cell cnv line number: $cnv_stem_path.");
            return 0;
        }
    }
}

1;
