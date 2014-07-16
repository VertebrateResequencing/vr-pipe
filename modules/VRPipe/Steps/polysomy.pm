
=head1 NAME

VRPipe::Steps::polysomy - a step

=head1 DESCRIPTION

This step runs polysomy algorithm which looks for acquired trisomy/tetrasomy in
a group of samples, then plots per-sample copy numbers for each chromosome.

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

class VRPipe::Steps::polysomy with VRPipe::StepRole {
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools exe',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_query_options => VRPipe::StepOption->create(
                description   => 'options to bcftools query, excluding -f and -s. bcftools is run to extract B-allele frequencies from input VCF files in order to be used by polysomy algorithm.',
                optional      => 1,
                default_value => "-r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
            ),
            control_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key to check on the input files to extract the sample identifier of the control from',
                default_value => 'sample_control'
            ),
            polysomy_exe => VRPipe::StepOption->create(
                description   => 'path to your polysomy exe',
                optional      => 1,
                default_value => 'polysomy'
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => '1 or more VCF files'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $cmk           = $options->{control_metadata_key};
            my $bcftools_exe  = $options->{bcftools_exe};
            my $bcftools_opts = $options->{bcftools_query_options};
            my $polysomy_exe  = $options->{polysomy_exe};
            
            if ($bcftools_opts =~ /\s-[fs]\s/) {
                $self->throw("bcftools_options should not include -f or -s");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => 0,
                    summary => "bcftools query merged.vcf.gz -f'[%CHROM\\t%BAF\\n]' -s \$control_sample $bcftools_opts \$vcf | polysomy -o \$outdir"
                )
            );
            
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                my @samples;
                my $fh = $vcf->openr;
                while (<$fh>) {
                    chomp;
                    next unless $_ =~ /^#CHROM/;
                    my @header = split(/\t/, $_);
                    @samples = @header[9 .. $#header];
                    last;
                }
                close($fh);
                
                my $meta = $vcf->metadata;
                
                my $control = $meta->{$cmk};
                $self->throw($vcf->path . " lacks a sample value for the $cmk metadata key") unless $control;
                
                my $req = $self->new_requirements(memory => 1000, time => 1);
                foreach my $query (@samples) {
                    my $sub_dir = $query;
                    
                    my $dist_file = $self->output_file(sub_dir => $sub_dir, output_key => 'dist_file', basename => 'dist.dat', type => 'txt', metadata => $meta);
                    $dist_file->add_metadata({ sample => $query });
                    
                    my $plot_file = $self->output_file(sub_dir => $sub_dir, output_key => 'plot_file', basename => 'dist.py', type => 'txt', metadata => $meta);
                    $plot_file->add_metadata({ sample => $query });
                    
                    my @outfiles  = ($plot_file, $dist_file);
                    my $vcf_path  = $vcf->path;
                    my $dist_path = $dist_file->path;
                    
                    my $this_cmd = "use VRPipe::Steps::polysomy; VRPipe::Steps::polysomy->run_and_check(vcf => q[$vcf_path], dist => q[$dist_path], bcftools => q[$bcftools_exe], bcftools_opts => q[$bcftools_opts], query => q[$query], polysomy => q[$polysomy_exe]);";
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
                
                }
            
            }
        };
    }
    
    method outputs_definition {
        return {
            dist_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => 1,
                description => 'output dist.dat file from polysomy',
            ),
            plot_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => 1,
                description => 'output dist.py file from polysomy',
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Detect trisomy/tetrasomy using Illumina's B-allele frequency (BAF) stored in a VCF file.";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method run_and_check (ClassName|Object $self: Str|File :$vcf!, Str|File :$dist!, Str :$bcftools!, Str :$bcftools_opts!, Str :$query!, Str :$polysomy!) {
        my $vcf_file  = VRPipe::File->get(path => $vcf);
        my $dist_file = VRPipe::File->get(path => $dist);
        
        my $cmd_line = "$bcftools query -f'[%CHROM\\t%BAF\\n]' -s $query $bcftools_opts " . $vcf_file->path . " | $polysomy -o " . $dist_file->dir;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $self->throw("File $dist doesn't exist..") unless -e $dist;
        
        my @chrs = `awk '{if(\$1==\"CN\" && \$3!=2.0)print \$2}' $dist`;
        if (@chrs) {
            chomp(@chrs);
            my $chrs_str = $query . ":" . join(",", @chrs);
            $vcf_file->merge_metadata({ polysomy_chrs => $chrs_str });
        }
        
        return 1;
    }

}

1;
