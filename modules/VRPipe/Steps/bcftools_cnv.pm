
=head1 NAME

VRPipe::Steps::bcftools_cnv - a step

=head1 DESCRIPTION

This step runs bcftools cnv on a VCF file which should always contain a 
control sample. For each sample in the VCF, the algorithm calls copy number 
variants compared to the control samples and generates the corresponding  copy
number difference plots by chromosome. If run after polysomy.pm, the  step
generates extra plots for chromosomes with evidence of trisomy/tetrasomy  not
picked up by bcftools cnv.

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

class VRPipe::Steps::bcftools_cnv with VRPipe::StepRole {
    method options_definition {
        return {
            bcftools_cnv_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools cnv exe',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_cnv_options => VRPipe::StepOption->create(
                description   => 'options to bcftools cnv, excluding -o, -c and -s',
                optional      => 1,
                default_value => '-p2 -S X,Y,MT'
            ),
            control_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key to check on the input files to extract the sample identifier of the control from',
                default_value => 'sample_control'
            ),
            python_exe => VRPipe::StepOption->create(
                description   => 'path to your python exe',
                optional      => 1,
                default_value => 'python'
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => '1 or more VCF files'
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $cmk              = $options->{control_metadata_key};
            my $bcftools_cnv_exe = $options->{bcftools_cnv_exe};
            my $bcftools_opts    = $options->{bcftools_cnv_options};
            my $python_exe       = $options->{python_exe};
            
            if ($bcftools_opts =~ /\s-[ocs]\s/) {
                $self->throw("bcftools_options should not include -o, -c or -s");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => 0,
                    summary => "bcftools cnv -c \$control_sample -s \$query_sample $bcftools_opts -o \$output_dir \$vcf_file"
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
                
                my $req = $self->new_requirements(memory => 1000, time => 5);
                
                foreach my $query (@samples) {
                    next if ($query eq $control);
                    
                    my $sub_dir = $query;
                    
                    my $summary_file = $self->output_file(sub_dir => $sub_dir, output_key => 'summary_file', basename => "summary.tab", type => 'txt', metadata => $meta);
                    $summary_file->add_metadata({ sample => $query });
                    
                    my $plot_file = $self->output_file(sub_dir => $sub_dir, output_key => 'plot_file', basename => "plot.$control.$query.py", type => 'txt', metadata => $meta);
                    
                    my @outfiles = ($summary_file, $plot_file);
                    
                    my @chroms = (1 .. 22, ("X", "Y", "MT"));
                    foreach my $chr (@chroms) {
                        my $png_file = $self->output_file(sub_dir => $sub_dir, output_key => 'png_files', basename => "plot.$control.$query.chr$chr.png", type => 'png', metadata => $meta);
                        push(@outfiles, $png_file);
                    }
                    
                    foreach my $sample (($control, $query)) {
                        $self->output_file(sub_dir => $sub_dir, basename => "summary.$sample.tab", type => 'txt', temporary => 1);
                        $self->output_file(sub_dir => $sub_dir, basename => "dat.$sample.tab",     type => 'txt', temporary => 1);
                        $self->output_file(sub_dir => $sub_dir, basename => "cn.$sample.tab",      type => 'txt', temporary => 1);
                    }
                    
                    my $vcf_path     = $vcf->path;
                    my $summary_path = $summary_file->path;
                    my $plot_path    = $plot_file->path;
                    
                    my $this_cmd = "use VRPipe::Steps::bcftools_cnv; VRPipe::Steps::bcftools_cnv->call_and_plot(vcf => q[$vcf_path], summary => q[$summary_path], plot => q[$plot_path], bcftools => q[$bcftools_cnv_exe], python => q[$python_exe], bcftools_opts => q[$bcftools_opts], control => q[$control], query => q[$query]);";
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
                
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            summary_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => 1,
                description => 'output summary file from bcftools cnv',
            ),
            plot_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => 1,
                description => 'output plot.py file from bcftools cnv',
            ),
            png_files => VRPipe::StepIODefinition->create(
                type            => 'png',
                min_files       => 0,
                max_files       => -1,
                description     => 'output plots by chr from bcftools cnv',
                check_existence => 0,
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Make copy number variation calls on a VCF using bcftools cnv.";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method call_and_plot (ClassName|Object $self: Str|File :$vcf!, Str|File :$summary!, Str|File :$plot!, Str :$bcftools!, Str :$python!, Str :$bcftools_opts!, Str :$control!, Str :$query!) {
        my $vcf_file     = VRPipe::File->get(path => $vcf);
        my $summary_file = VRPipe::File->get(path => $summary);
        my $plot_file    = VRPipe::File->get(path => $plot);
        
        my $cmd_line = "$bcftools cnv -c $control -s $query -o " . $summary_file->dir . " $bcftools_opts " . $vcf_file->path;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my @chrs = ();
        foreach ($vcf_file->meta_value('polysomy_chrs')) {
            my $chrs_str;
            if ($_ =~ /$query:/) {
                ($chrs_str) = $_ =~ m/$query:(.*)/;
            }
            elsif ($_ =~ /$control:/) {
                ($chrs_str) = $_ =~ m/$control:(.*)/;
            }
            my @cs = split(/,/, $chrs_str);
            push(@chrs, @cs);
        }
        
        if (@chrs) {
            my %seen = ();
            my @unique_chrs = grep { !$seen{$_}++ } @chrs;
            foreach my $chr (@unique_chrs) {
                my $cmd_line = "$python " . $plot_file->path . " -c $chr";
                system($cmd_line) && $self->throw("failed to run [$cmd_line]");
            }
        }
        
        return 1;
    }

}

1;
