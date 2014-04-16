
=head1 NAME

VRPipe::Steps::bcftools_generate_sites_file - a step

=head1 DESCRIPTION

Generates sites files from vcf.gz or bcf, in same dir as the input file.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bcftools_generate_sites_file extends VRPipe::Steps::bcftools {
    around options_definition {
        return {
            %{ $self->$orig },
            'genotypes_file' => VRPipe::StepOption->create(
                description => 'absolute path to vcf or bcf from which to generate sites file; optional if you supply this via the datasource input',
                optional    => 1,
            )
        };
    }
    
    method inputs_definition {
        return {
            genotypes_bcf => VRPipe::StepIODefinition->create(
                type        => 'bcf',
                min_files   => 0,
                max_files   => 1,
                description => 'bcf file containing calls at all desired sites'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            
            my $genotypes_file_path;
            if ($options->{genotypes_file}) {
                $genotypes_file_path = file($options->{genotypes_file});
                $self->throw("genotypes_file must be an absolute path") unless $genotypes_file_path->is_absolute;
            }
            else {
                my ($genotypes_bcf) = @{ $self->inputs->{genotypes_bcf} || [] };
                $genotypes_file_path = $genotypes_bcf->path if $genotypes_bcf;
            }
            $genotypes_file_path || $self->throw("A genotypes file is required, either via the genotypes_file option or as an input from the datasource");
            
            my $basename = $genotypes_file_path->basename;
            $basename =~ s/\.[bv]cf(.gz)?$/.tab.txt/;
            my $sites_file = $self->output_file(
                output_key => 'sites_file',
                output_dir => $genotypes_file_path->dir->stringify, # goes into the vcf dir
                basename   => $basename, type => 'txt'
            );
            
            my $output_path = $sites_file->path;
            
            my $cmd = "$bcftools_exe view $genotypes_file_path | grep -v ^# | cut -f1,2 > $output_path";
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch_wrapped_cmd('VRPipe::Steps::bcftools_generate_sites_file', 'generate_sites', [$cmd, $req, { output_files => [$sites_file], block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return { sites_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'genomic sites file', max_files => 1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generates sites files from vcf.gz or bcf.";
    }
    
    method max_simultaneous {
        return 0;                                                   # meaning unlimited
    }
    
    method generate_sites (ClassName|Object $self: Str $cmd_line) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($input_path, $output_path) = $cmd_line =~ /view (\S+) .* > (\S+)$/;
        
        my $input_file = VRPipe::File->create(path => $input_path);
        my $input_recs = $input_file->num_records;
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_recs = $output_file->num_records;
        
        unless ($output_recs == $input_recs) {
            $output_file->unlink;
            $self->throw("Output sites file has different number of data lines from input vcf(input $input_recs, output $output_recs)");
        }
        else {
            return 1;
        }
    }
}

1;
