
=head1 NAME

VRPipe::Steps::vcf_combine_genotypes_and_sites - a step

=head1 DESCRIPTION

Merge a sites only file with the genotypes found in another file,  keeping the
FILTER, ID and INFO fields found in the sites only file.

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::vcf_combine_genotypes_and_sites with VRPipe::StepRole {
    method options_definition {
        return {
            sites_vcf_path => VRPipe::StepOption->create(description => 'Absolute path to sites vcf'),
        };
    }
    
    method inputs_definition {
        return {
            genotype_vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', description => 'One or more bgzipped and indexed vcf files containing the genotype data which will be used in the output vcf', max_files => -1, metadata => { seq_no => 'a sequence number assigned by the split for reassembly in correct order' }),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $options = $self->options;
            
            my $sites_vcf_path = $options->{sites_vcf_path};
            
            # create temporary fofn of files to merge
            my $concat_list      = $self->output_file(basename => "concat_list.txt", type => 'txt', temporary => 1);
            my $concat_list_path = $concat_list->path;
            my @sorted_vcf_files = sort { $a->metadata->{seq_no} <=> $b->metadata->{seq_no} } @{ $self->inputs->{genotype_vcf_files} };
            $concat_list->create_fofn(\@sorted_vcf_files);
            
            # define output file
            my $combined_meta     = $self->common_metadata($self->inputs->{genotype_vcf_files});
            my $combined_vcf      = $self->output_file(output_key => 'combined_vcf', basename => "combined.vcf.gz", type => 'vcf', metadata => $combined_meta);
            my $combined_vcf_path = $combined_vcf->path;
            
            my $chrom = exists $$combined_meta{chrom} ? qq[, chrom => q[$$combined_meta{chrom}]] : '';
            my $cmd = "use VRPipe::Steps::vcf_combine_genotypes_and_sites; VRPipe::Steps::vcf_combine_genotypes_and_sites->combine_genotypes_and_sites(output_vcf => q[$combined_vcf_path], sites_vcf => q[$sites_vcf_path], concat_list => q[$concat_list_path]$chrom);";
            $self->dispatch_vrpipecode($cmd, $req, { output_files => [$combined_vcf] });
        };
    }
    
    method outputs_definition {
        return { combined_vcf => VRPipe::StepIODefinition->create(type => 'vcf', description => 'a merged vcf file', max_files => 1) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merge a sites only VCF file with the genotypes found in another VCF file, keeping the FILTER, ID and INFO fields found in the sites only file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    # 2012-09-11: This method is modified for VRPipe from Petr Danecek's run-mpileup
    method combine_genotypes_and_sites (ClassName|Object $self: Str|File :$output_vcf!, Str|File :$sites_vcf!, Str|File :$concat_list!, Str :$chrom = '.') {
        my $output_file = VRPipe::File->get(path => file($output_vcf));
        my $concat_file = VRPipe::File->get(path => file($concat_list));
        my $sites_file = VRPipe::File->create(path => file($sites_vcf));
        
        $output_file->disconnect;
        
        open(my $sites_fh, "tabix -h $sites_vcf $chrom |") or $self->throw("tabix -h $sites_vcf $chrom: $!");
        open(my $out_fh,   "| bgzip -c > $output_vcf")     or $self->throw("bgzip -c > $output_vcf: $!");
        my $files_fh = $concat_file->openr;
        my ($gt_fh, $gt_file, $gt_line, $header_printed, @columns, @header, $column_map);
        while (my $site_line = <$sites_fh>) {
            if (substr($site_line, 0, 1) eq '#') {
                if ($header_printed) { next; }
                if (substr($site_line, 1, 1) eq '#') { push @header, $site_line; } # this header contains the filters
                next;
            }
            while (!defined $gt_line) {
                if (defined $gt_fh) {
                    close($gt_fh) or $self->throw("close tabix -h $gt_file $chrom");
                    $gt_fh = undef;
                }
                $gt_file = <$files_fh>;
                if (!defined $gt_file) { last; }
                chomp($gt_file);
                open($gt_fh, "tabix -h $gt_file $chrom |") or $self->throw("tabix -h $gt_file $chrom: $!");
                while ($gt_line = <$gt_fh>) {
                    if (substr($gt_line, 0, 1) ne '#') { last; } # data line
                    if (substr($gt_line, 1, 1) eq '#') { next; } # header line
                    
                    # column line
                    my @cols = split(/\t/, $gt_line);
                    if (!$header_printed) {
                        @columns = @cols;
                        print $out_fh @header;
                        print $out_fh $gt_line;
                        $header_printed = 1;
                        next;
                    }
                    $column_map = $self->pad_columns(\@columns, \@cols);
                }
            }
            if (!defined $gt_line) {
                $self->throw("Out of sync: no gt_line in $gt_file after $site_line");
            }
            my $site_idx = 0;
            my $gt_idx   = 0;
            for (my $i = 0; $i < 5; $i++) {
                $site_idx = index($site_line, "\t", $site_idx + 1);
            }
            if (substr($site_line, 0, $site_idx) ne substr($gt_line, 0, $site_idx)) {
                # Either out of sync or an indel was modified
                my @site_items = split(/\t/, $site_line);
                my @site_alts  = split(/,/,  $site_items[4]);
                for (my $i = 0; $i < 5; $i++) { $gt_idx = index($gt_line, "\t", $gt_idx + 1); }
                my @gt_items = split(/\t/, substr($gt_line, 0, $gt_idx));
                my @gt_alts = split(/,/, $gt_items[4]);
                if (@site_alts != @gt_alts) { $self->throw("Out of sync: $site_idx\n$site_line$gt_line"); }
                my $site_ref_len = length($site_items[3]);
                my $gt_ref_len   = length($gt_items[3]);
                if ($site_ref_len >= $gt_ref_len) { $self->throw("Out of sync: $site_idx\n$site_line$gt_line"); }
                my $diff = substr($gt_items[3], $site_ref_len);
                
                for (my $i = 0; $i < @site_alts; $i++) {
                    if ($gt_alts[$i] ne $site_alts[$i] . $diff) { $self->throw("Out of sync: $site_idx, [$site_alts[$i]] [$gt_alts[$i].$diff]\n$site_line$gt_line"); }
                }
            }
            else {
                $gt_idx = $site_idx;
            }
            print $out_fh substr($site_line, 0, -1);
            for (my $i = 0; $i < 3; $i++) {
                $gt_idx = index($gt_line, "\t", $gt_idx + 1);
            }
            if (!defined $column_map) {
                print $out_fh substr($gt_line, $gt_idx);
            }
            else {
                my @cols = split(/\t/, substr($gt_line, $gt_idx + 1));
                chomp($cols[-1]);
                for (my $i = 8; $i < @$column_map; $i++) {
                    my $idx = $$column_map[$i];
                    if   ($idx == -1) { print $out_fh "\t."; }
                    else              { print $out_fh "\t", $cols[$idx - 8]; }
                }
                print $out_fh "\n";
            }
            $gt_line = <$gt_fh>;
        }
        $concat_file->close;
        close($gt_fh)    or $self->throw("close tabix -h $gt_file $chrom");
        close($out_fh)   or $self->throw("close bgzip -c > $output_vcf");
        close($sites_fh) or $self->throw("close tabix -h $sites_vcf $chrom");
        
        $output_file->update_stats_from_disc;
        my $output_records = $output_file->num_records;
        
        my $input_records = 0;
        $files_fh = $concat_file->openr;
        while (my $input_path = <$files_fh>) {
            chomp $input_path;
            my $input_file = VRPipe::File->get(path => $input_path);
            $input_records += $input_file->num_records;
        }
        $concat_file->close;
        
        unless ($output_records == $input_records) {
            $output_file->unlink;
            $self->throw("Output VCF has different number of data lines from inputs (inputs $input_records, output $output_records)");
        }
        else {
            return 1;
        }
    }
    
    method pad_columns (ClassName|Object $self: ArrayRef $cols1!, ArrayRef $cols2!) {
        if (@$cols1 < @$cols2) { $self->throw(sprintf "Not ready for this, sorry, expected fewer columns (%d!<%d)", @$cols1, @$cols2); }
        my $has1 = {};
        my $has2 = {};
        for (my $i = 0; $i < @$cols1; $i++) { $$has1{ $$cols1[$i] } = $i; }
        my $need_map = 0;
        for (my $i = 0; $i < @$cols2; $i++) {
            if (!exists($$has1{ $$cols2[$i] })) { $self->throw("The column [$$cols2[$i]] not seen previously."); }
            if ($i != $$has1{ $$cols2[$i] }) { $need_map = 1; }
            $$has2{ $$cols2[$i] } = $i;
        }
        if (@$cols1 == @$cols2 && !$need_map) { return undef; }
        my @map;
        for (my $i = 0; $i < @$cols1; $i++) {
            my $cname = $$cols1[$i];
            push @map, exists($$has2{$cname}) ? $$has2{$cname} : -1;
        }
        return \@map;
    }
}

1;
