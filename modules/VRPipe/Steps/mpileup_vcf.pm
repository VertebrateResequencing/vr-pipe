
=head1 NAME

VRPipe::Steps::mpileup_vcf - a step

=head1 DESCRIPTION

Runs samtools mpileup and bcftools for one or more BAMS, generating one VCF
without an intermediate BCF

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

class VRPipe::Steps::mpileup_vcf extends VRPipe::Steps::mpileup_bcf {
    around options_definition {
        return { %{ $self->$orig },
                 samtools_mpileup_options => VRPipe::StepOption->create(description => 'samtools mpileup options excluding -f and -b. Since this will be piped into bcftools view, it is recommended that the -u option is set.',                                                                              optional => 1, default_value => '-DSV -C50 -m2 -F0.0005 -d 10000 -ug'),
                 bcftools_exe             => VRPipe::StepOption->create(description => 'path to bcftools executable',                                                                                                                                                                                          optional => 1, default_value => 'bcftools'),
                 bcftools_view_options    => VRPipe::StepOption->create(description => 'bcftools view options',                                                                                                                                                                                                optional => 1, default_value => '-p 0.99 -vcgN'),
                 sample_sex_file          => VRPipe::StepOption->create(description => 'File listing the sex (M or F) of samples. If not provided, will call on all samples in the bcf header. If provided, calls will be made on the intersection of the samples in the file and samples in the bcf header.', optional => 1),
                 assumed_sex              => VRPipe::StepOption->create(description => 'If M or F is not present for a sample in the sample sex file, then this sex is assumed',                                                                                                                               optional => 1, default_value => 'F') };
    }
    
    method body_sub {
        return sub {
            my $self            = shift;
            my $options         = $self->options;
            my $samtools        = $options->{samtools_exe};
            my $reference_fasta = $options->{reference_fasta};
            my $mpileup_opts    = $options->{samtools_mpileup_options};
            my $bcftools        = $options->{bcftools_exe};
            my $bcf_view_opts   = $options->{bcftools_view_options};
            my $assumed_sex     = $options->{assumed_sex};
            my $sample_sex_file;
            
            if ($options->{sample_sex_file}) {
                $sample_sex_file = Path::Class::File->new($options->{sample_sex_file});
                $self->throw("sample_sex_file must be an absolute path") unless $sample_sex_file->is_absolute;
            }
            if ($self->inputs->{sites_file}) {
                $self->throw("bcftools_view_options cannot contain the -l option if a sites_file is an input to this step") if ($bcf_view_opts =~ /-l/);
                my $sites_file = $self->inputs->{sites_file}[0];
                $bcf_view_opts .= " -l " . $sites_file->path;
            }
            
            my $bams_list = $self->output_file(basename => "bams.list", type => 'txt', temporary => 1);
            my $bams_list_path = $bams_list->path;
            $bams_list->create_fofn($self->inputs->{bam_files});
            my $vcf_meta = $self->common_metadata($self->inputs->{bam_files});
            my @bam_ids = map { $_->id } @{ $self->inputs->{bam_files} };
            $vcf_meta->{caller} = 'samtools_mpileup_bcftools';
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my $female_ploidy = 2;
            my $male_ploidy   = 2;
            
            my $element_meta = $self->step_state->dataelement->result;
            my $basename;
            if (defined $element_meta->{chrom}) {
                my $chrom = $element_meta->{chrom};
                my $from  = $element_meta->{from};
                my $to    = $element_meta->{to};
                $female_ploidy = $element_meta->{female_ploidy} if defined $element_meta->{female_ploidy};
                $male_ploidy   = $element_meta->{male_ploidy}   if defined $element_meta->{male_ploidy};
                my $region = "${chrom}_${from}-${to}";
                $$vcf_meta{chrom}  = $chrom;
                $$vcf_meta{from}   = $from;
                $$vcf_meta{to}     = $to;
                $$vcf_meta{seq_no} = $element_meta->{chunk_id};
                my $override_file = $element_meta->{chunk_override_file};
                my $override      = do $override_file;
                
                if (exists $override->{"$region"}->{samtools_exe}) {
                    $samtools = $override->{"$region"}->{samtools_exe};
                }
                if (exists $override->{"$region"}->{samtools_mpileup_options}) {
                    $mpileup_opts = $override->{"$region"}->{samtools_mpileup_options};
                }
                if (exists $override->{"$region"}->{bcftools_exe}) {
                    $bcftools = $override->{"$region"}->{bcftools_exe};
                }
                if (exists $override->{"$region"}->{bcftools_view_options}) {
                    $bcf_view_opts = $override->{"$region"}->{bcftools_view_options};
                }
                $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe     => 'samtools',
                                                                      version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                                                                      summary => "samtools mpileup $mpileup_opts -r \$region -f \$reference_fasta \$bam_files | bcftools view $bcf_view_opts -s \$samples_file - | bgzip -c > \$vcf_file"));
                $basename = "$region.mpileup";
                $mpileup_opts .= " -r $chrom:$from-$to";
            }
            else {
                $basename = 'mpileup';
                $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe     => 'samtools',
                                                                      version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                                                                      summary => "samtools mpileup $mpileup_opts -f \$reference_fasta \$bam_files | bcftools view $bcf_view_opts -s \$samples_file - | bgzip -c > \$vcf_file"));
            }
            
            my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => $basename . '.vcf.gz', type => 'vcf', metadata => $vcf_meta);
            my $temp_samples_path = $self->output_file(basename => $basename . '.samples', type => 'txt', temporary => 1)->path;
            
            my $mpileup_cmd  = qq[$samtools mpileup $mpileup_opts -f $reference_fasta -b $bams_list_path];
            my $bcftools_cmd = qq[$bcftools view $bcf_view_opts -s $temp_samples_path];
            my $cmd_line     = qq[$mpileup_cmd | $bcftools_cmd - | bgzip -c > ] . $vcf_file->path;
            
            my $args = qq['$cmd_line', '$temp_samples_path', source_file_ids => [qw(@bam_ids)], female_ploidy => '$female_ploidy', male_ploidy => '$male_ploidy', assumed_sex => '$assumed_sex'];
            $args .= qq[, sample_sex_file => '$sample_sex_file'] if $sample_sex_file;
            my $cmd = "use VRPipe::Steps::bcf_to_vcf; VRPipe::Steps::bcf_to_vcf->bcftools_call_with_sample_file($args);";
            $self->dispatch_vrpipecode($cmd, $req, { output_files => [$vcf_file] });
        };
    }
    
    method description {
        return "Run samtools mpileup and bcftools for one or more bams, generating one vcf without an intermediate bcf";
    }
    
    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a vcf file for each set of one or more input bams') };
    }
    
    # method run_mpileup (ClassName|Object $self: Str $cmd_line, Int $min_recs) {
    #     system($cmd_line) && $self->throw("failed to run [$cmd_line]");
    #
    #     my ($output_path) = $cmd_line =~ /> (\S+)$/;
    #     my $output_file = VRPipe::File->get(path => $output_path);
    #     $output_file->update_stats_from_disc;
    #     my $output_recs = $output_file->num_records;
    #
    #     if ($output_recs < $min_recs) {
    #         $output_file->unlink;
    #         $self->throw("Output VCF has $output_recs data lines, fewer than expected minimum $min_recs");
    #     } else {
    #         return 1;
    #     }
    # }
}

1;
