
=head1 NAME

VRPipe::Steps::vcf_create_merged_sites_list - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::vcf_create_merged_sites_list with VRPipe::StepRole {
    method options_definition {
        return { vcf_isec_exe => VRPipe::StepOption->create(description => 'path to your vcf-isec executable', optional => 1, default_value => 'vcf-isec'), };
    }
    
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', description => 'compressed vcf files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options  = $self->options;
            my $isec_exe = $options->{vcf_isec_exe};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my @input_set;
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                push @input_set, $vcf_file->path;
            }
            
            my $site_list = $self->output_file(output_key => 'merged_sites_list', basename => 'merged.list', type => 'txt');
            
            my $output_path = $site_list->path;
            my $this_cmd    = "$isec_exe -f -n +1 @input_set | awk '!(\\\$1~/^#/)' | cut -f 1,2 | uniq > $output_path";
            
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_create_merged_sites_list', 'merge_sites_list', [$this_cmd, $req, { output_files => [$site_list] }]);
        };
    }
    
    method outputs_definition {
        return { merged_sites_list => VRPipe::StepIODefinition->create(type => 'txt', description => 'site list containing', max_files => 1, check_existence => 0), };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges compressed VCFs, creating a single output VCF";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method merge_sites_list (ClassName|Object $self: Str $cmd_line) {
        my ($first_input_path, $output_path) = $cmd_line =~ /^\S+ -f -n \+1 (\S+) .* (\S+)$/;
        
        my $first_input_file = VRPipe::File->get(path => $first_input_path);
        my $first_input_records = $first_input_file->num_records;
        
        $first_input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        if ($output_lines < $first_input_records) {
            $output_file->unlink;
            $self->throw("Output sites list has $output_lines, fewer than number of records in first input $first_input_records");
        }
        else {
            return 1;
        }
    }
}

1;
