
=head1 NAME

VRPipe::Steps::vcf_merge_different_samples - a step

=head1 DESCRIPTION

Run bcftools merge to merge a set of VCFs each with a different sample,
creating a single output VCF with all the input samples.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2014 Genome Research Limited.

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

class VRPipe::Steps::vcf_merge_different_samples with VRPipe::StepRole {
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools executable',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_options => VRPipe::StepOption->create(
                description   => 'Options for the bcftools merge command. Does not support the --print-header option',
                optional      => 1,
                default_value => '-O z'
            )
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
            
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            my $bcfopts      = $options->{bcftools_options};
            
            my @input_set;
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                push @input_set, $vcf_file->path;
            }
            my $merged_basename = 'merged.vcf.gz';
            
            my $merged_meta = $self->common_metadata($self->inputs->{vcf_files});
            my $merged_vcf  = $self->output_file(output_key => 'merged_vcf', basename => $merged_basename, type => 'vcf', metadata => $merged_meta);
            my $output_path = $merged_vcf->path;
            
            if (@input_set == 1) {
                # merge doesn't work on 1 input file; just symlink the input to
                # output
                my $source = $self->inputs->{vcf_files}->[0];
                $source->symlink($merged_vcf);
                return 1;
            }
            
            my $this_cmd = qq[$bcftools_exe merge $bcfopts @input_set > $output_path];
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => VRPipe::StepCmdSummary->determine_version($bcftools_exe, '^Version: (.+)$'),
                    summary => "bcftools merge $bcfopts \@input_vcfs > \$output_path"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_merge_different_samples', 'merge_vcf', [$this_cmd, $req, { output_files => [$merged_vcf] }]);
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
        return "Merges compressed VCFs using bcftools merge which contain different samples to produce a single VCF containing all the input samples";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method merge_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($first_input_path, $output_path) = $cmd_line =~ /^.+ (\S+?\.gz) .* (\S+)$/;
        
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
