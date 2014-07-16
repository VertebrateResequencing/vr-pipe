
=head1 NAME

VRPipe::Steps::plot_polysomy - a step

=head1 DESCRIPTION

This step runs plot-polysomy to generate per-sample copy numbers for each 
chromosome.

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

class VRPipe::Steps::plot_polysomy with VRPipe::StepRole {
    method options_definition {
        return {
            plot_polysomy_exe => VRPipe::StepOption->create(
                description   => 'path to your plot-polysomy script',
                optional      => 1,
                default_value => 'plot-polysomy.py'
            ),
        };
    }
    
    method inputs_definition {
        return {
            dist_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => '1 or more dist.dat files from polysomy',
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $plot_polysomy_exe = $options->{plot_polysomy_exe};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'plot-polysomy',
                    version => 0,
                    summary => "plot-polysomy -o \$png_file \@[title:dir]"
                )
            );
            
            my @args;
            foreach my $file (@{ $self->inputs->{dist_files} }) {
                my $dist  = $file->metadata;
                my $query = $dist->{sample};
                $query =~ s/:/_/g; #create safe string
                push(@args, "$query:" . $file->dir);
            }
            
            my $vcf = $self->inputs->{vcf_files};
            
            my $merged_meta = $self->combined_metadata($self->inputs->{dist_files});
            my $png_file    = $self->output_file(output_key => 'png_file', basename => 'copy_numbers.png', type => 'png', metadata => $merged_meta);
            my $req         = $self->new_requirements(memory => 1000, time => 1);
            my $cmd_line    = "$plot_polysomy_exe -o " . $png_file->path . " @args";
            $self->dispatch([$cmd_line, $req]);
        };
    }
    
    method outputs_definition {
        return {
            png_file => VRPipe::StepIODefinition->create(
                type        => 'png',
                min_files   => 1,
                max_files   => 1,
                description => 'output copy number plot from polysomy',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run plot-polysomy to generate copy number plots by chromosome.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}

1;
