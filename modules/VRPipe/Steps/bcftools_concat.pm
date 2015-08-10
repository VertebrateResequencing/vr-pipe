
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
                default_value => ''
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
                description => 'after bcftools concat, option to pipe output vcf through another command, e.g. "vcf-annotate --fill-ICF" to fill AC, AN, and ICF annotations',
                optional    => 1
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'var',
                max_files   => -1,
                description => 'vcf files to concat',
                metadata    => {
                    chrom    => 'chromosome',
                    from     => 'region start',
                    to       => 'region end',
                    optional => ['chrom', 'from', 'to']
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
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
            
            my %seen;
            my @chroms = map { $_->metadata->{chrom} ? $_->metadata->{chrom} : () } @{ $self->inputs->{vcf_files} };
            my @uniq_chroms = grep { !$seen{$_}++ } @chroms;
            if ($merge_by_chr && @uniq_chroms) {
                foreach my $chr (@uniq_chroms) {
                    my @vcf_files = ();
                    foreach my $vcf_file (@{ $self->inputs->{vcf_files} }) {
                        if ($vcf_file->metadata->{chrom} eq $chr) {
                            push(@vcf_files, $vcf_file);
                        }
                    } # create temporary fofn of files to merge
                    my $merge_list = $self->output_file(basename => "merge_list.chr$chr.txt", type => 'txt', temporary => 1);
                    my @sorted_vcf_files = sort { $a->metadata->{from} cmp $b->metadata->{from} } @vcf_files;
                    my $file_list_id    = VRPipe::FileList->create(files => \@sorted_vcf_files)->id;
                    my $concat_meta     = $self->common_metadata(\@vcf_files);
                    my $concat_vcf      = $self->output_file(output_key => 'concat_vcf', basename => "merged.chr$chr.vcf.gz", type => 'vcf', metadata => $concat_meta);
                    my $vcf_index       = $self->output_file(output_key => 'vcf_index', basename => "merged.chr$chr.vcf.gz.tbi", type => 'bin', metadata => $concat_meta);
                    my $merge_list_path = $merge_list->path;
                    my $concat_vcf_path = $concat_vcf->path;
                    my $cmd             = qq[($bcftools_exe concat$opts -f $merge_list_path$cut$filter | bgzip -c > $concat_vcf_path) && $bcftools_exe index -ft $concat_vcf_path];
                    my $this_cmd        = "use VRPipe::Steps::bcftools_concat; VRPipe::Steps::bcftools_concat->concat_and_check(q[$cmd], input_file_list => $file_list_id);";
                    $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$concat_vcf, $merge_list, $vcf_index] });
                }
            }
            else {
                my $merge_list = $self->output_file(basename => "merge_list.txt", type => 'txt', temporary => 1);
                
                my $file_list_id;
                if (@chroms) {
                    my @sorted_vcf_files = sort { $a->metadata->{chrom} cmp $b->metadata->{chrom} || $a->metadata->{from} <=> $b->metadata->{from} } @{ $self->inputs->{vcf_files} };
                    $file_list_id = VRPipe::FileList->create(files => \@sorted_vcf_files)->id;
                }
                else {
                    $file_list_id = VRPipe::FileList->create(files => $self->inputs->{vcf_files})->id;
                }
                my $concat_meta     = $self->common_metadata($self->inputs->{vcf_files});
                my $concat_vcf      = $self->output_file(output_key => 'concat_vcf', basename => "merged.vcf.gz", type => 'vcf', metadata => $concat_meta);
                my $vcf_index       = $self->output_file(output_key => 'vcf_index', basename => "merged.vcf.gz.tbi", type => 'bin', metadata => $concat_meta);
                my $merge_list_path = $merge_list->path;
                my $concat_vcf_path = $concat_vcf->path;
                my $cmd             = qq[($bcftools_exe concat$opts -f $merge_list_path$cut$filter | bgzip -c > $concat_vcf_path) && $bcftools_exe index -ft $concat_vcf_path];
                my $this_cmd        = "use VRPipe::Steps::bcftools_concat; VRPipe::Steps::bcftools_concat->concat_and_check(q[$cmd], input_file_list => $file_list_id);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$concat_vcf, $merge_list, $vcf_index] });
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
    
    method concat_and_check (ClassName|Object $self: Str $cmd_line!, Int :$input_file_list!) {
        my ($input_path)  = $cmd_line =~ /-f (\S+)/;
        my ($output_path) = $cmd_line =~ /(\S+)$/;
        $input_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $output_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $out_file = VRPipe::File->get(path => $output_path);
        
        my @input_files = VRPipe::FileList->get(id => $input_file_list)->files;
        my $fofn = VRPipe::File->get(path => $input_path);
        $fofn->create_fofn(\@input_files);
        
        $out_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
    }

}

1;
