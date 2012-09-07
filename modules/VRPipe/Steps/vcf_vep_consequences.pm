
=head1 NAME

VRPipe::Steps::vcf_vep_consequences - a step

=head1 DESCRIPTION

Runs the vcf2consequences_vep program, which applies Consequence annotations
generated using the Ensembl VEP to VCF files, and optionally Grantham Matrix
and GERP scores.

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

class VRPipe::Steps::vcf_vep_consequences with VRPipe::StepRole {
    method options_definition {
        return {
            'vcf2consequences_options' => VRPipe::StepOption->create(description => 'options to vcf2consequences_vep, excluding -v and -i'),
            'vcf2consequences_exe'     => VRPipe::StepOption->create(
                description   => 'path to your vcf2consequences executable',
                optional      => 1,
                default_value => 'vcf2consequences_vep'
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'annotated vcf files',
                max_files   => -1
            ),
            vep_txt => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'vep analysis output file',
                metadata    => { source_vcf => 'the vcf file analysed by the VEP' },
                max_files   => -1
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options  = $self->options;
            my $con_exe  = $options->{'vcf2consequences_exe'};
            my $con_opts = $options->{'vcf2consequences_options'};
            
            if ($con_opts =~ /-[v,i]/) {
                $self->throw("vcf2consequences_options should not include the -i or -v option");
            }
            
            my $req = $self->new_requirements(memory => 5000, time => 1);
            
            # put vep file metadata into hash for vcf name lookup
            my (%vep_files);
            foreach my $vep (@{ $self->inputs->{vep_txt} }) {
                my $path       = $vep->path->stringify;
                my $vep_meta   = $vep->metadata;
                my $source_vcf = $vep_meta->{source_vcf};
                $vep_files{$source_vcf} = $path;
            }
            
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                my $input_path = $vcf_file->path;
                my $vep_txt_path = $vep_files{$input_path} || $self->throw("got no recal file for $input_path");
                
                my $basename = $vcf_file->basename;
                $basename =~ s/\.vcf(.gz)?$/.conseq.vcf.gz/;
                my $conseq_vcf = $self->output_file(output_key => 'conseq_vcf', basename => $basename, type => 'vcf', metadata => $vcf_file->metadata);
                
                my $output_path = $conseq_vcf->path;
                
                my $this_cmd = "$con_exe -v $input_path -i $vep_txt_path $con_opts | bgzip -c > $output_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_vep_consequences', 'consequence_vcf', [$this_cmd, $req, { output_files => [$conseq_vcf] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            conseq_vcf => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'annotated vcf file with VEP consequences',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generates VCF files annotated with VEP consequences";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method consequence_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $output_path) = $cmd_line =~ /\S+ -v (\S+) .* bgzip -c > (\S+)$/;
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
