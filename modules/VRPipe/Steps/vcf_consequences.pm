
=head1 NAME

VRPipe::Steps::vcf_consequences - a step

=head1 DESCRIPTION

Runs the old vcf2consequences program aganinst one or more VCF files to generate VCFs annotated with consequence, Superceded by vep_analysis / vcf2consequences_vep, which uses the Ensembl VEP and is much quicker.[Should probably get rid of this Step]

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

class VRPipe::Steps::vcf_consequences with VRPipe::StepRole {
    method options_definition {
        return { 'vcf2consequences_options' => VRPipe::StepOption->create(description => 'options to vcf2consequences, excluding -v'),
                 'vcf2consequences_exe'     => VRPipe::StepOption->create(
                                                                      description   => 'path to your vcf2consequences executable',
                                                                      optional      => 1,
                                                                      default_value => 'vcf2consequences'),
                 'tabix_exe' => VRPipe::StepOption->create(description   => 'path to your tabix executable',
                                                           optional      => 1,
                                                           default_value => 'tabix') };
    }
    
    method inputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type        => 'vcf',
                                                               description => 'annotated vcf files',
                                                               max_files   => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options   = $self->options;
            my $tabix_exe = $options->{tabix_exe};
            my $con_exe   = $options->{'vcf2consequences_exe'};
            my $con_opts  = $options->{'vcf2consequences_options'};
            
            if ($con_opts =~ /-v/) {
                $self->throw("vcf2consequences_options should not include the reference or -v option");
            }
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                my $basename = $vcf_file->basename;
                if ($basename =~ /\.vcf.gz$/) {
                    $basename =~ s/\.vcf.gz$/.conseq.vcf.gz/;
                }
                else {
                    $basename =~ s/\.vcf$/.conseq.vcf/;
                }
                my $conseq_vcf = $self->output_file(output_key => 'conseq_vcf', basename => $basename,          type => 'vcf');
                my $tbi        = $self->output_file(output_key => 'tbi_file',   basename => $basename . '.tbi', type => 'bin');
                
                my $input_path  = $vcf_file->path;
                my $output_path = $conseq_vcf->path;
                
                my $this_cmd = "$con_exe -v $input_path $con_opts | bgzip -c > $output_path; $tabix_exe -f -p vcf $output_path;";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_consequences', 'consequence_vcf', [$this_cmd, $req, { output_files => [$conseq_vcf, $tbi] }]);
            }
        };
    }
    
    method outputs_definition {
        return { conseq_vcf => VRPipe::StepIODefinition->create(type        => 'vcf',
                                                                description => 'annotated vcf file with consequences',
                                                                max_files   => -1),
                 tbi_file => VRPipe::StepIODefinition->create(type        => 'bin',
                                                              description => 'a tbi file',
                                                              max_files   => -1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Adds consequence annotation to VCF files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method consequence_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /\S+ -v (\S+) .* vcf (\S+);$/;
        my $input_file = VRPipe::File->get(path => $input_path);
        
        my $input_lines = $input_file->lines;
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        # Should have a few extra header lines
        unless ($output_lines >= $input_lines) {
            $output_file->unlink;
            $self->throw("Output VCF has $output_lines lines, less than input $input_lines");
        }
        else {
            return 1;
        }
    }
}

1;
