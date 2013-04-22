
=head1 NAME

VRPipe::Steps::vcf_sites - a step

=head1 DESCRIPTION

Generates sites files from VCF, in same dir as VCF.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::vcf_sites with VRPipe::StepRole {
    method options_definition {
        return { 'genotypes_vcf' => VRPipe::StepOption->create(description => 'absolute path to vcf from which to generate sites file') };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options   = $self->options;
            my $vcf_path = $options->{genotypes_vcf};
            
            my $vcf = Path::Class::File->new($vcf_path);
            my $basename = $vcf->basename;

            my $cat_exe = $basename =~ /\.vcf.gz$/ ? 'zcat' : 'cat';

            $basename =~ s/\.vcf(.gz)?$/.tab.txt/;
            my $sites_file = $self->output_file(output_key => 'sites_file',
                output_dir => $vcf->dir->stringify,     # goes into the vcf dir
                basename => $basename, type => 'txt');

            my $output_path = $sites_file->path;
            
            my $cmd = "$cat_exe $vcf_path | grep -v ^# | cut -f1,2 > $output_path";
            my $req = $self->new_requirements(memory => 500, time => 1);
            $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_sites', 'gen_vcf_sites', [$cmd, $req, { output_files => [$sites_file] }]);
        };
    }
    
    method outputs_definition {
        return { sites_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'genomic sites file', max_files => 1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generates sites files from VCF";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

    method gen_vcf_sites (ClassName|Object $self: Str $cmd_line) {

        system($cmd_line) && $self->throw("failed to run [$cmd_line]");

        my ($cat_exe, $input_path, $output_path) = $cmd_line =~ /^(\S+) (\S+) .* > (\S+)$/;

        my $input_recs=0;
        my $pipe = "$cat_exe $input_path |";
        my $fh;
        open($fh, $pipe) || $self->throw("Couldn't open '$pipe': $!");
        while (<$fh>) {
            next if /^#/;
            $input_recs++;
        }
        close($fh);
        
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
