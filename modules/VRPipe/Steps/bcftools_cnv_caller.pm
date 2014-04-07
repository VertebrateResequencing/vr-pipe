
=head1 NAME

VRPipe::Steps::bcftools_cnv_caller - a step

=head1 DESCRIPTION

This step runs bcftools_cnv_caller, which calls loss of heterozygosity on
genotype array data in VCF format.

*** more documentation to come

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

class VRPipe::Steps::bcftools_cnv_caller with VRPipe::StepRole {
    method options_definition {
        return {
            bcftools_cnv_caller_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools_cnv_caller script',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_cnv_caller_options => VRPipe::StepOption->create(
                description   => 'options to bcftools_cnv_caller, excluding -o, -c and -s',
                optional      => 1,
                default_value => '-p0'
            ),
            control_metadata_key => VRPipe::StepOption->create(
                description   => 'the metadata key to check on the input files to extract the sample identifier of the control from',
                default_value => 'control'
            ),
            cnv_analysis_type => VRPipe::StepOption->create(
                description   => 'type of cnv analysis, added to file metadata for downstream processing',
                optional      => 1,
                default_value => 'bcftools_cnv'
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
            
            my $cmk               = $options->{control_metadata_key};
            my $caller            = $options->{bcftools_cnv_caller_exe};
            my $caller_opts       = $options->{bcftools_cnv_caller_options};
            my $cnv_analysis_type = $options->{cnv_analysis_type};
            if ($caller_opts =~ /\s-[ocs]\s/) {
                $self->throw("bcftools_cnv_caller_options should not include -o, -c or -s");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools_cnv_caller',
                    version => 0,
                    summary => "bcftools_cnv_caller cnv -c \$control_sample -s \$query_sample $caller_opts -o \$output_dir \$vcf_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 2000, time => 1);
            
            foreach my $vcf (@{ $self->inputs->{vcf_files} }) {
                print "vcf path:", $vcf->path . "\n";
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
                
                my $new_meta = { cnv_analysis_type => $cnv_analysis_type };
                $vcf->add_metadata($new_meta);
                
                my $meta = $vcf->metadata;
                
                print "one: @samples\n";
                
                my $control = $meta->{$cmk};
                $self->throw($vcf->path . " lacks a sample value for the $cmk metadata key") unless $control;
                
                foreach my $query (@samples) {
                    next if ($query eq $control);
                    
                    my $sub_dir = "$control-$query";
                    
                    print "two\t$sub_dir\n";
                    
                    my $summary_file = $self->output_file(sub_dir => $sub_dir, output_key => 'bcftools_cnv_txt_files', basename => "summary.tab", type => 'txt', metadata => $meta);
                    
                    my @chroms = (1 .. 22, ("X", "Y", "MT"));
                    foreach my $chr (@chroms) {
                        $self->output_file(sub_dir => $sub_dir, output_key => 'bcftools_cnv_png_files', basename => "cnv-$query-$control-chr$chr.png", type => 'png', metadata => $meta);
                    }
                    
                    $self->output_file(sub_dir => $sub_dir, output_key => 'bcftools_cnv_txt_files', basename => 'baf.tab',  type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, output_key => 'bcftools_cnv_txt_files', basename => 'dat.tab',  type => 'txt', temporary => 1);
                    $self->output_file(sub_dir => $sub_dir, output_key => 'bcftools_cnv_txt_files', basename => 'plots.py', type => 'txt', temporary => 1);
                    
                    my $this_cmd = "$caller cnv -c $control -s $query -o " . $summary_file->dir . " $caller_opts " . $vcf->path;
                    
                    $self->dispatch([$this_cmd, $req])
                
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            bcftools_cnv_png_files => VRPipe::StepIODefinition->create(
                type            => 'png',
                min_files       => 0,
                max_files       => 30,
                description     => 'output from the bcftools cnv',
                check_existence => 0
            ),
            bcftools_cnv_txt_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => 1,
                description => 'output from the bcftools cnv',
            )
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

}

1;
