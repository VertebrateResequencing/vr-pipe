
=head1 NAME

VRPipe::Steps::vcf_multi_filter - a step

=head1 DESCRIPTION

Soft-filters input VCFs in parallel from a multi-column datasource, using multiple filter programs and option files for each VCF. 
The filter programs and their accociated option files are each provided as parameters delimited by '#'. 

vcf-filter#vcf-annotate
vcf_filter_opt_file#vcf_annotate_opt_file

There needs to be one program and filter options file for each column in the datasource.

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::vcf_multi_filter extends VRPipe::Steps::vcf_filter {
    method options_definition {
        return { 'vcf-filter_programs' => VRPipe::StepOption->create(description => 'path to one or more filter executables, one per vcf to be filtered in parallel, delimited by #'),
                 'vcf-filter_files'    => VRPipe::StepOption->create(description => 'one or more filter program option files, one per filter executable, delimited by #'),
                 'tabix_exe'           => VRPipe::StepOption->create(description => 'path to your tabix executable', optional => 1, default_value => 'tabix'), };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options         = $self->options;
            my @filter_programs = split('#', $options->{'vcf-filter_programs'});
            my @filter_files    = split('#', $options->{'vcf-filter_files'});
            my $tabix_exe       = $options->{tabix_exe};
            my $cat_exe;
            
            my $idx = 0;
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                my $basename = $vcf_file->basename;
                if ($basename =~ /\.vcf.gz$/) {
                    $basename =~ s/\.vcf.gz$/.filt.vcf.gz/;
                    $cat_exe = 'zcat';
                }
                else {
                    $basename =~ s/\.vcf$/.filt.vcf/;
                    $cat_exe = 'cat';
                }
                my $filter_exe  = $filter_programs[$idx];
                my $filter_file = $filter_files[$idx];
                $idx++;
                
                my $filtered_vcf = $self->output_file(output_key => 'filtered_vcf', basename => $basename,          type => 'vcf');
                my $tbi          = $self->output_file(output_key => 'tbi_file',     basename => $basename . '.tbi', type => 'bin');
                
                my $input_path  = $vcf_file->path;
                my $output_path = $filtered_vcf->path;
                
                my $this_cmd = "$cat_exe $input_path | $filter_exe -f $filter_file | bgzip -c > $output_path; $tabix_exe -f -p vcf $output_path";
                
                $self->dispatch_wrapped_cmd('VRPipe::Steps::vcf_filter', 'filter_vcf', [$this_cmd, $req, { output_files => [$filtered_vcf, $tbi] }]);
            }
        };
    }
    
    method description {
        return "Soft-filters input VCFs in parallel from a multi-column datasource, using multiple filter programs and option files";
    }
}

1;
