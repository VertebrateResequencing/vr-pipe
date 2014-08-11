
=head1 NAME

VRPipe::Steps::create_imputation_ref_panel - a step

=head1 DESCRIPTION

This step converts reference VCF files into haps/legend/sample formats which 
are used for imputation reference panels by impute2 and shapeit. It takes bed 
files as input which determine the coordinates of the chunks of haplotypes to 
be generated.

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

class VRPipe::Steps::create_imputation_ref_panel with VRPipe::StepRole {
    method options_definition {
        return {
            vcf_gensample_exe => VRPipe::StepOption->create(
                description   => 'path to your vcf-gensample exe',
                optional      => 1,
                default_value => 'vcf-gensample'
            ),
            ref_vcf => VRPipe::StepOption->create(
                description => 'path to a reference VCF file which may contain the string "{CHROM}" that will be expanded to chromomome names (more than one VCF can be specified each of which should be placed in a separate directory where VRPipe must have write permission).',
                optional    => 1,
            ),
            min_AC => VRPipe::StepOption->create(
                description   => 'skip sites with the number of alternate alleles smaller than min_AC',
                optional      => 0,
                default_value => 0,
            ),
        };
    }
    
    method inputs_definition {
        return {
            bed_files => VRPipe::StepIODefinition->create(
                type        => 'txt',
                min_files   => 1,
                max_files   => -1,
                description => '1 or more bed files'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $vcf_gensample_exe = $options->{vcf_gensample_exe};
            my $ref_vcf           = $options->{ref_vcf};
            my $min_AC            = $options->{min_AC};
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'vcf-gensample',
                    version => 0,
                    summary => "vcf-gensample -m \$min_ac -r \$chr:\$from-\$to -o outdir/\$from-\$to \$vcf_file",
                )
            );
            
            if ($ref_vcf) {
                my %dirs = ();
                my @ref_vcfs = split(/\s+/, $ref_vcf);
                foreach my $str (@ref_vcfs) {
                    my ($d, $f) = split(/\/([^\/]+)$/, $str);
                    if (exists $dirs{$d}) {
                        self->throw("Reference VCF files must be placed in separate directories");
                    }
                    else {
                        $dirs{$d} = ".";
                    }
                }
                
                foreach my $bed_file (@{ $self->inputs->{bed_files} }) {
                    my $fh = $bed_file->openr;
                    while (<$fh>) {
                        chomp;
                        next if ($_ =~ /^#/);
                        my @chunk = split(/\t/, $_);
                        my $chr   = $chunk[0];
                        my $from  = $chunk[1];
                        my $to    = $chunk[2];
                        foreach my $ref (@ref_vcfs) {
                            my $ref_path = $self->expand_chrom($ref, $chr);
                            my $ref_file = VRPipe::File->create(path => $ref_path);
                            my $outdir = $ref_file->dir . "/$chr";
                            
                            my $haps_file   = $self->output_file(output_dir => "$outdir", output_key => 'haplotype_files', basename => "$from-${to}_haps.gz",   type => 'bin', metadata => { chrom => $chr, from => $from, to => $to });
                            my $legend_file = $self->output_file(output_dir => "$outdir", output_key => 'haplotype_files', basename => "$from-${to}_legend.gz", type => 'bin', metadata => { chrom => $chr, from => $from, to => $to });
                            my $sample_file = $self->output_file(output_dir => "$outdir", output_key => 'haplotype_files', basename => "$from-${to}_samples",   type => 'txt', metadata => { chrom => $chr, from => $from, to => $to });
                            my $sites_file  = $self->output_file(output_dir => "$outdir", output_key => 'sites_files',     basename => "$from-$to.ref.sites",   type => 'txt');
                            
                            my @outfiles = ($haps_file, $legend_file, $sample_file, $sites_file);
                            
                            my $req = $self->new_requirements(memory => 1000, time => 1);
                            my $this_cmd = "use VRPipe::Steps::create_imputation_ref_panel; VRPipe::Steps::create_imputation_ref_panel->create_chunks(chunk => [qw(@chunk)], ref_vcf => q[$ref_path], vcf_gensample => q[$vcf_gensample_exe], outdir => q[$outdir], min_AC => q[$min_AC]);";
                            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles, block_and_skip_if_ok => 1 });
                        }
                    
                    }
                    $fh->close;
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            haplotype_files => VRPipe::StepIODefinition->create(
                type            => 'any',
                min_files       => 0,
                max_files       => -1,
                description     => 'reference haplotype chunks in haps/legend/sample format',
                check_existence => 0
            ),
            sites_files => VRPipe::StepIODefinition->create(
                type            => 'txt',
                min_files       => 0,
                max_files       => -1,
                description     => 'reference haplotype positions in chunks',
                check_existence => 0
            ),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Split and convert reference VCF files into smaller haps/legend/sample chunks";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method expand_chrom (Str $path, Str $region) {
        if (!defined $region) { return $path; }
        $region =~ s/:.*$//;
        $path =~ s/{CHROM}/$region/g;
        return $path;
    }
    
    method create_chunks (ClassName|Object $self: Str|File :$ref_vcf!, ArrayRef[Str] :$chunk!, Str :$vcf_gensample!, Str :$outdir!, Str :$min_AC!){
        my @region = @{$chunk};
        my $chr    = $region[0];
        my $from   = $region[1];
        my $to     = $region[2];
        
        my $cmd_line = "$vcf_gensample -m $min_AC -r $chr:$from-$to $ref_vcf -l -o $outdir/$from-$to";
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        # additionally print out the sites in the reference panel
        my $legend_file = VRPipe::File->get(path => "$outdir/$from-${to}_legend.gz");
        my $sites_file  = VRPipe::File->get(path => "$outdir/$from-$to.ref.sites");
        my $in          = $legend_file->openr;
        my $out         = $sites_file->openw;
        while (my $line = <$in>) {
            my @items = split(/\s+/, $line);
            print $out $items[1], "\n";
        }
        $out->close;
        $in->close;
        
        return 1;
    }

}

1;
