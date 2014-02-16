
=head1 NAME

VRPipe::Steps::vcf_merge - a step

=head1 DESCRIPTION

Run vcf-isec to merge a set of  VCFs, creating a single output VCF

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::vcf_merge with VRPipe::StepRole {
    method options_definition {
        return {
            'vcf-isec_exe'      => VRPipe::StepOption->create(description => 'path to your vcf-isec executable',                                                                                                                                                                                         optional      => 1,          default_value => 'vcf-isec'),
            vcf_isec_options    => VRPipe::StepOption->create(description => 'Options for the vcf-isec command. Does not support the --prefix (-p) option',                                                                                                                                              default_value => '-f -n +1', optional      => 1),
            metadata_priority   => VRPipe::StepOption->create(description => 'String of the form key#value1:value2:value3 etc where each of the input vcf files has the "key" metadata. Input vcfs will be ordered according to the values listed. Values must be unique amongst the step input files.', optional      => 1),
            post_merge_vcftools => VRPipe::StepOption->create(description => 'After merging with vcf-isec, option to pipe output vcf through a vcftools command, e.g. "vcf-annotate --fill-ICF" to fill AC, AN, and ICF annotations',                                                                    optional      => 1),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'compressed vcf files',
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options           = $self->options;
            my $isec_exe          = $options->{'vcf-isec_exe'};
            my $isec_options      = $options->{vcf_isec_options};
            my $metadata_priority = $options->{metadata_priority};
            my $post_filter       = $options->{post_merge_vcftools};
            
            my ($key, $priority_list, @values);
            if ($metadata_priority) {
                ($key, $priority_list) = split(/#/, $metadata_priority);
                @values = split(/:/, $priority_list);
            }
            
            my @ordered_vcfs;
            my $merged_basename = '';
            if ($metadata_priority) {
                my %vcf_paths;
                foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                    my $value = $vcf->metadata->{$key};
                    $self->throw("More than one input vcf has the metadata $value for key $key") if (exists $vcf_paths{$value});
                    $vcf_paths{$value} = $vcf;
                }
                foreach my $value (@values) {
                    $self->throw("None of the input vcfs have the metadata $value for key $key") unless (exists $vcf_paths{$value});
                    push @ordered_vcfs, $vcf_paths{$value};
                }
            }
            else {
                @ordered_vcfs = @{ $self->inputs->{vcf_files} };
            }
            
            my @input_set;
            foreach my $vcf_file (@ordered_vcfs) {
                push @input_set, $vcf_file->path;
                my $basename = $vcf_file->basename;
                $basename =~ s/vcf.gz$//;
                $merged_basename .= $basename;
            }
            $merged_basename .= 'merged.vcf.gz';
            
            my $merged_meta = $self->common_metadata($self->inputs->{vcf_files});
            my $merged_vcf = $self->output_file(output_key => 'merged_vcf', basename => $merged_basename, type => 'vcf', metadata => $merged_meta);
            
            my $filter = $post_filter ? " | $post_filter" : '';
            
            my $output_path = $merged_vcf->path;
            my $this_cmd    = qq[$isec_exe $isec_options @input_set$filter | bgzip -c > $output_path];
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'vcf-isec',
                    version => 0,
                    summary => "vcf-isec $isec_options \@input_vcfs$filter | bgzip -c > \$output_path"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_merge', 'merge_vcf', [$this_cmd, $req, { output_files => [$merged_vcf] }]);
        };
    }
    
    method outputs_definition {
        return {
            merged_vcf => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'a merged vcf file',
                max_files   => 1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges compressed VCFs using vcf-isec which contain the same set of samples, creating a single output VCF with that set of samples";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method merge_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($first_input_path, $output_path) = $cmd_line =~ /^\S+.*? (\S+?.gz) .* (\S+)$/;
        
        my $first_input_file = VRPipe::File->get(path => $first_input_path);
        my $first_input_lines = $first_input_file->lines;
        
        $first_input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        if ($output_lines < $first_input_lines) {
            $output_file->unlink;
            $self->throw("Output VCF has $output_lines, fewer than first input $first_input_lines");
        }
        else {
            return 1;
        }
    }
}

1;
