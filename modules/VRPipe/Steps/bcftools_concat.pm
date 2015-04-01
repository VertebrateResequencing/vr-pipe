
=head1 NAME

VRPipe::Steps::bcftools_concat - a step

=head1 DESCRIPTION

Runs bcftools concat against an input set of VCFs with 'chrom', 'from' and 'to'
metadata needed for sorting purposes. The step combines the input VCF files 
generating a single merged VCF or one merged VCF per chromosome.

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

class VRPipe::Steps::bcftools_concat with VRPipe::StepRole {
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to bcftools executable',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_concat_opts => VRPipe::StepOption->create(
                description   => 'options to bcftools concat (excluding -f, -o and -O)',
                optional      => 1,
                default_value => '--ligate'
            ),
            concat_sites_only => VRPipe::StepOption->create(
                description   => 'do not output genotype information to the concatenated vcf file',
                optional      => 1,
                default_value => 0
            ),
            merge_by_chromosome => VRPipe::StepOption->create(
                description   => 'concatenate chunks belonging to each chromosome separately, generating one merged VCF for each chromosome',
                optional      => 1,
                default_value => 0
            ),
            post_concat_opts => VRPipe::StepOption->create(
                description => 'after vcf-concat, option to pipe output vcf through a vcftools command, e.g. "vcf-annotate --fill-ICF" to fill AC, AN, and ICF annotations',
                optional    => 1
            ),
            tabix_exe => VRPipe::StepOption->create(
                description   => 'path to your tabix exe',
                optional      => 1,
                default_value => 'tabix'
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                max_files   => -1,
                description => 'vcf files to concat',
                metadata    => {
                    chrom => 'chromosome',
                    from  => 'region start',
                    to    => 'region end'
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $tabix_exe     = $options->{tabix_exe};
            my $bcftools_exe  = $options->{bcftools_exe};
            my $bcftools_opts = $options->{bcftools_concat_opts};
            my $sites_only    = $options->{concat_sites_only};
            my $merge_by_chr  = $options->{merge_by_chromosome};
            my $post_filter   = $options->{post_concat_opts};
            
            if ($bcftools_opts =~ /\s-[fOo]\s/) {
                $self->throw("bcftools_options should not include -f, -o or -O");
            }
            
            my $opts   = $bcftools_opts ? " $bcftools_opts" : '';
            my $cut    = $sites_only    ? ' | cut -f 1-8'   : '';
            my $filter = $post_filter   ? " | $post_filter" : '';
            my $req = $self->new_requirements(memory => 5000, time => 1);
            
            if ($merge_by_chr) {
                my %seen;
                my @meta_chrs = map  { $_->metadata->{chrom} } @{ $self->inputs->{vcf_files} };
                my @chrs      = grep { !$seen{$_}++ } @meta_chrs;
                foreach my $chr (@chrs) {
                    my @vcf_files = ();
                    foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                        if ($vcf_file->metadata->{chrom} eq $chr) {
                            push(@vcf_files, $vcf_file);
                        }
                    } # create temporary fofn of files to merge
                    my $merge_list = $self->output_file(basename => "merge_list.chr$chr.txt", type => 'txt', temporary => 1);
                    my @sorted_vcf_files = sort { $a->metadata->{from} <=> $b->metadata->{from} } @vcf_files;
                    $merge_list->create_fofn(\@sorted_vcf_files);
                    my $concat_meta     = $self->common_metadata(\@vcf_files);
                    my $concat_vcf      = $self->output_file(output_key => 'concat_vcf', basename => "merged.chr$chr.vcf.gz", type => 'vcf', metadata => $concat_meta);
                    my $vcf_index       = $self->output_file(output_key => 'vcf_index', basename => "merged.chr$chr.vcf.gz.tbi", type => 'bin');
                    my $merge_list_path = $merge_list->path;
                    my $concat_vcf_path = $concat_vcf->path;
                    my $cmd             = qq[($bcftools_exe concat$opts -f $merge_list_path$cut$filter | bgzip -c > $concat_vcf_path) && $tabix_exe -f -p vcf $concat_vcf_path];
                    $self->dispatch([$cmd, $req, { output_files => [$concat_vcf, $merge_list, $vcf_index] }]);
                }
            }
            else {
                my $merge_list = $self->output_file(basename => "merge_list.txt", type => 'txt', temporary => 1);
                my @sorted_vcf_files = sort { $a->metadata->{chrom} <=> $b->metadata->{chrom} || $a->metadata->{from} <=> $b->metadata->{from} } @{ $self->inputs->{vcf_files} };
                $merge_list->create_fofn(\@sorted_vcf_files);
                my $concat_meta     = $self->common_metadata($self->inputs->{vcf_files});
                my $concat_vcf      = $self->output_file(output_key => 'concat_vcf', basename => "merged.vcf.gz", type => 'vcf', metadata => $concat_meta);
                my $vcf_index       = $self->output_file(output_key => 'vcf_index', basename => "merged.vcf.gz.tbi", type => 'bin');
                my $merge_list_path = $merge_list->path;
                my $concat_vcf_path = $concat_vcf->path;
                my $cmd             = qq[($bcftools_exe concat$opts -f $merge_list_path$cut$filter | bgzip -c > $concat_vcf_path) && $tabix_exe -f -p vcf $concat_vcf_path];
                $self->dispatch([$cmd, $req, { output_files => [$concat_vcf, $merge_list, $vcf_index] }]);
            }
        
        };
    }
    
    method outputs_definition {
        return {
            concat_vcf => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                min_files   => 1,
                max_files   => -1,
                description => 'Concatenated VCF files'
            ),
            vcf_index => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'index of the concatanated vcf file',
                min_files   => 1,
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run bcftools to combine or concatenate input VCF chunks";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }

}

1;
