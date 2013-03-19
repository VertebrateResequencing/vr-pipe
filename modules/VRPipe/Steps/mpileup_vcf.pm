
=head1 NAME

VRPipe::Steps::mpileup_vcf - a step

=head1 DESCRIPTION

Runs samtools mpileup and bcftools for one or more BAMS, generating one VCF
without an intermediate BCF

=head1 AUTHOR

Chris Joyce    <cj5@sanger.ac.uk>. Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::mpileup_vcf extends VRPipe::Steps::bcf_to_vcf {
    around options_definition {
        return {
            %{ $self->$orig },
            samtools_exe             => VRPipe::StepOption->create(description => 'path to samtools executable',                                                                                                             optional => 1, default_value => 'samtools'),
            samtools_mpileup_options => VRPipe::StepOption->create(description => 'samtools mpileup options excluding -f and -b. Since this will be piped into bcftools view, it is recommended that the -u option is set.', optional => 1, default_value => '-DSV -C50 -m2 -F0.0005 -d 10000 -ug'),
            reference_fasta          => VRPipe::StepOption->create(description => 'absolute path to reference genome fasta')
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(type => 'bam', max_files => -1, description => '1 or more bam files to call variants'),
            bai_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => -1,
                description => 'index files for the input bam files'
            ),
            sites_file => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0, max_files => 1, description => 'Optional sites file for calling only at the given sites'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{bam_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($vcf_meta);
            
            my $samtools        = $options->{samtools_exe};
            my $bcftools        = $options->{bcftools_exe};
            my $mpileup_opts    = $options->{samtools_mpileup_options};
            my $bcf_view_opts   = $options->{bcftools_view_options};
            my $assumed_sex     = $options->{assumed_sex};
            my $minimum_records = $options->{minimum_records};
            my $post_filter     = $options->{post_calling_vcftools};
            
            my $reference_fasta = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $reference_fasta->is_absolute;
            
            if ($mpileup_opts =~ /-f|-b|$reference_fasta/) {
                $self->throw("samtools_mpileup_options should not include the reference, the -f or the -b options");
            }
            
            my $sample_sex_file;
            if ($options->{sample_sex_file}) {
                $sample_sex_file = file($options->{sample_sex_file});
                $self->throw("sample_sex_file must be an absolute path") unless $sample_sex_file->is_absolute;
            }
            if ($self->inputs->{sites_file}) {
                $self->throw("bcftools_view_options cannot contain the -l option if a sites_file is an input to this step") if ($bcf_view_opts =~ /-l/);
                my $sites_file = $self->inputs->{sites_file}[0];
                $bcf_view_opts .= " -l " . $sites_file->path;
            }
            my $filter = $post_filter ? " | $post_filter" : '';
            
            my $bams_list_path = $self->output_file(basename => "bams.list", type => 'txt', temporary => 1)->path;
            my @bam_ids = map { $_->id } @{ $self->inputs->{bam_files} };
            $vcf_meta->{caller} = 'samtools_mpileup_bcftools';
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my $summary_opts = $mpileup_opts;
            my $basename     = 'mpileup';
            my ($female_ploidy, $male_ploidy) = (2, 2);
            if (defined $$vcf_meta{chrom} && defined $$vcf_meta{from} && defined $$vcf_meta{to}) {
                my ($chrom, $from, $to) = ($$vcf_meta{chrom}, $$vcf_meta{from}, $$vcf_meta{to});
                $summary_opts .= ' -r $region';
                $mpileup_opts .= " -r $chrom:$from-$to";
                $basename      = "${chrom}_${from}-${to}.$basename";
                $female_ploidy = $vcf_meta->{female_ploidy} if defined $vcf_meta->{female_ploidy};
                $male_ploidy   = $vcf_meta->{male_ploidy} if defined $vcf_meta->{male_ploidy};
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'samtools',
                    version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                    summary => "samtools mpileup $summary_opts -f \$reference_fasta -b \$bams_list | bcftools view $bcf_view_opts -s \$samples_file -$filter | bgzip -c > \$vcf_file"
                )
            );
            
            my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => $basename . '.vcf.gz', type => 'vcf', metadata => $vcf_meta);
            my $temp_samples_path = $self->output_file(basename => $basename . '.samples', type => 'txt', temporary => 1)->path;
            
            my $mpileup_cmd  = qq[$samtools mpileup $mpileup_opts -f $reference_fasta -b $bams_list_path];
            my $bcftools_cmd = qq[$bcftools view $bcf_view_opts -s $temp_samples_path];
            my $cmd_line     = qq[$mpileup_cmd | $bcftools_cmd -$filter | bgzip -c > ] . $vcf_file->path;
            
            my $args = qq['$cmd_line', '$temp_samples_path', source_file_ids => [qw(@bam_ids)], female_ploidy => '$female_ploidy', male_ploidy => '$male_ploidy', assumed_sex => '$assumed_sex'];
            $args .= qq[, sample_sex_file => '$sample_sex_file'] if $sample_sex_file;
            my $cmd = "use VRPipe::Steps::bcf_to_vcf; VRPipe::Steps::bcf_to_vcf->bcftools_call_with_sample_file($args, minimum_records => $minimum_records);";
            $self->dispatch_vrpipecode($cmd, $req, { output_files => [$vcf_file] });
        };
    }
    
    method description {
        return "Run samtools mpileup and bcftools for one or more bams, generating one vcf without an intermediate bcf";
    }
    
    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a vcf file for each set of one or more input bams') };
    }
}

1;
