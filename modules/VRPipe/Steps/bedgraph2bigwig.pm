
=head1 NAME

VRPipe::Steps::bedgraph2bigwig - a step

=head1 DESCRIPTION

The step converts input bedGraph files to bigWig format.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::bedgraph2bigwig with VRPipe::StepRole {
    method options_definition {
        return {
            bedGraphToBigWig_exe => VRPipe::StepOption->create(
                description   => 'path to your bedGraphToBigWig binary',
                optional      => 1,
                default_value => 'bedGraphToBigWig',
            ),
            bedGraphToBigWig_options => VRPipe::StepOption->create(
                description => 'options to bedGraphToBigWig',
                optional    => 1,
            ),
            reference_index => VRPipe::StepOption->create(
                description => 'absolute path to the fasta index (.fai) file associated with reference fasta file',
            ),
        };
    }
    
    method inputs_definition {
        return {
            bdg_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'bedGraph files',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self             = shift;
            my $options          = $self->options;
            my $bedGraphToBigWig = $options->{bedGraphToBigWig_exe};
            my $fai_file         = file($options->{reference_index});
            $self->throw("reference_index must be an absolute path") unless $fai_file->is_absolute;
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => "$bedGraphToBigWig",
                    version => '',
                    summary => "$bedGraphToBigWig \$input_bedGraph chrom.sizes \$output_BigWig"
                )
            );
            
            foreach my $input_bedGraph (@{ $self->inputs->{bdg_files} }) {
                my $output_BigWig = $self->output_file(output_key => 'bigwig_files', basename => $input_bedGraph->basename . ".bw", type => 'bin');
                my @outfiles      = ($output_BigWig);
                my $this_cmd      = qq[$bedGraphToBigWig ] . $input_bedGraph->path . qq[ $options->{reference_index} ] . $output_BigWig->path . qq [ $options->{bedGraphToBigWig_options}];
                my $req           = $self->new_requirements(memory => 2000, time => 1);
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bedgraph2bigwig', 'convert_and_check', [$this_cmd, $req, { output_files => \@outfiles }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bigwig_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => -1,
                description => 'BigWig files',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Converts bedGraph files to bigWig format";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method convert_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in) = $cmd_line =~ /bedGraphToBigWig (\S+)/;
        my $input_bedGraph = VRPipe::File->get(path => file($in)->absolute);
        
        $input_bedGraph->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
    }

}

1;
