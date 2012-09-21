
=head1 NAME

VRPipe::Steps::vcf_annotate - a step

=head1 DESCRIPTION

Generates annotated VCF files from one or more input VCFs using vcf-annotate. 
Optionally, one can make two passes with vcf-annotate.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::vcf_annotate with VRPipe::StepRole {
    method options_definition {
        return {
            'vcf-annotate_options'   => VRPipe::StepOption->create(description => 'vcf-annotate pass 1 options'),
            'vcf-annotate_2_options' => VRPipe::StepOption->create(description => 'vcf-annotate pass 2 options', optional => 1),
            'vcf-annotate_exe' => VRPipe::StepOption->create(
                description   => 'path to your vcf-annotate executable',
                optional      => 1,
                default_value => 'vcf-annotate'
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'vcf files',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options   = $self->options;
            my $an_exe    = $options->{'vcf-annotate_exe'};
            my $an_opts   = $options->{'vcf-annotate_options'};
            my $an_2_opts = $options->{'vcf-annotate_2_options'};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                my $basename = $vcf_file->basename;
                my $cat_exe = $basename =~ /\.vcf.gz$/ ? 'zcat' : 'cat';
                $basename =~ s/\.vcf(.gz)?$/.annot.vcf.gz/;
                
                my $annotated_vcf = $self->output_file(output_key => 'annotated_vcf', basename => $basename, type => 'vcf', metadata => $vcf_file->metadata);
                
                my $input_path  = $vcf_file->path;
                my $output_path = $annotated_vcf->path;
                
                my $annotate = "$an_exe $an_opts";
                if ($an_2_opts) {
                    $annotate .= " | $an_exe $an_2_opts";
                }
                my $this_cmd = "$cat_exe $input_path | $annotate | bgzip -c > $output_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_annotate', 'annotate_vcf', [$this_cmd, $req, { output_files => [$annotated_vcf] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            annotated_vcf => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'annotated vcf file',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Annotate VCF files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method annotate_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /^\S+ (\S+) .* bgzip -c > (\S+)$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        
        my $input_recs = $input_file->num_records;
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_recs = $output_file->num_records;
        
        unless ($output_recs == $input_recs) {
            $output_file->unlink;
            $self->throw("Output VCF has different number of data lines from input (input $input_recs, output $output_recs)");
        }
        else {
            return 1;
        }
    }
}

1;
