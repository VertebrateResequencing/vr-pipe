
=head1 NAME

VRPipe::Steps::combine_bcftools_cnvs - a step

=head1 DESCRIPTION

This step runs cmp-cnvs.pl to combine the summary results from bcftools cnv 
across the samples.

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::combine_bcftools_cnvs with VRPipe::StepRole {
    method options_definition {
        return {
            compare_cnv_script => VRPipe::StepOption->create(
                description   => 'path to your cmp-cnvs.pl which combines bcftools outputs - skips if not provided',
                optional      => 1,
                default_value => 'cmp-cnvs.pl',
            ),
            compare_cnv_options => VRPipe::StepOption->create(
                description   => 'options to cmp-cnvs.pl which combines output from bcftools',
                optional      => 1,
                default_value => '-q 2'
            ),
        };
    }
    
    method inputs_definition {
        return {
            summary_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'summary files from bcftools cnv',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $cmp_cnvs_exe  = $options->{compare_cnv_script};
            my $cmp_cnvs_opts = $options->{compare_cnv_options};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'cmp-cnvs.pl',
                    version => 0,
                    summary => "cmp-cnvs.pl $cmp_cnvs_opts \@dirs > combined_cnvs.txt"
                )
            );
            
            my @dirs;
            foreach my $file (@{ $self->inputs->{summary_files} }) {
                push(@dirs, $file->dir);
            }
            my $merged_meta = $self->combined_metadata($self->inputs->{summary_files});
            my $cmp_file = $self->output_file(output_key => 'cmp_file', basename => 'combined_cnvs.txt', type => 'txt', metadata => $merged_meta);
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            my $cmd_line = "$cmp_cnvs_exe $cmp_cnvs_opts @dirs > " . $cmp_file->path;
            $self->dispatch([$cmd_line, $req]);
        };
    }
    
    method outputs_definition {
        return {
            cmp_file => VRPipe::StepIODefinition->create(
                type            => 'txt',
                min_files       => 1,
                max_files       => 1,
                description     => 'combined results from cmp-cnvs.pl',
                check_existence => 0,
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Combine bcftools cnv summary results across the samples.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
