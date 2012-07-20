
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
        return { 'vcf-isec_exe' => VRPipe::StepOption->create(description   => 'path to your vcf-isec executable',
                                                              optional      => 1,
                                                              default_value => 'vcf-isec'),
                 'tabix_exe' => VRPipe::StepOption->create(description   => 'path to your tabix executable',
                                                           optional      => 1,
                                                           default_value => 'tabix') };
    }
    
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type        => 'vcf',
                                                               description => 'compressed vcf files',
                                                               max_files   => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options   = $self->options;
            my $tabix_exe = $options->{tabix_exe};
            my $isec_exe  = $options->{'vcf-isec_exe'};
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my @input_set;
            my $merged_basename = '';
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                push @input_set, $vcf_file->path;
                
                my $basename = $vcf_file->basename;
                $basename =~ s/vcf.gz$//;
                $merged_basename .= $basename;
            }
            $merged_basename .= 'merged.vcf.gz';
            
            my $merged_vcf = $self->output_file(output_key => 'merged_vcf', basename => $merged_basename,       type => 'vcf');
            my $tbi        = $self->output_file(output_key => 'tbi_file',   basename => "$merged_basename.tbi", type => 'bin');
            
            my $output_path = $merged_vcf->path;
            my $this_cmd    = "$isec_exe -f -n +1 @input_set | bgzip -c > $output_path; $tabix_exe -f -p vcf $output_path";
            
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_merge', 'merge_vcf', [$this_cmd, $req, { output_files => [$merged_vcf, $tbi] }]);
        };
    }
    
    method outputs_definition {
        return { merged_vcf => VRPipe::StepIODefinition->create(type        => 'vcf',
                                                                description => 'a merged vcf file',
                                                                max_files   => 1),
                 tbi_file => VRPipe::StepIODefinition->create(type        => 'bin',
                                                              description => 'a tbi file',
                                                              max_files   => -1) };
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
    
    method merge_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($first_input_path, $output_path) = $cmd_line =~ /^\S+ -f -n \+1 (\S+) .* (\S[^;]+);/;
        
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
