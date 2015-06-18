
=head1 NAME

VRPipe::Steps::polysomy - a step

=head1 DESCRIPTION

This step runs polysomy algorithm which looks for acquired trisomy/tetrasomy in
a group of samples, then plots per-sample copy numbers for each chromosome.

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, 2015 Genome Research Limited.

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
    use VRPipe::Schema;
    
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to your bcftools exe',
                optional      => 1,
                default_value => 'bcftools'
            ),
            bcftools_polysomy_options => VRPipe::StepOption->create(
                description   => 'options to bcftools polysomy, excluding -o and -s',
                optional      => 1,
                default_value => '-t ^MT,Y'
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
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $cmk           = $options->{control_metadata_key};
            my $bcftools_exe  = $options->{bcftools_exe};
            my $bcftools_opts = $options->{bcftools_polysomy_options};
            my $python_exe    = $options->{python_exe};
            
            if ($bcftools_opts =~ /\s-[os]\s/) {
                $self->throw("bcftools_polysomy_options should not include -o or -s");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => 0,
                    summary => "bcftools polysomy \$vcf -s \$control_sample $bcftools_opts -o \$outdir"
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
                my $sample_count;
                foreach my $query (@samples) {
                    # if a sample name is duplicated, it will have a number
                    # prefix in the vcf header
                    $sample_count++;
                    my $real_sample_name = $query;
                    $real_sample_name =~ s/^$sample_count://;
                    
                    my $sub_dir = $query;
                    
                    my $single_sample_meta = {};
                    while (my ($key, $val) = each %$meta) {
                        next unless defined $val;
                        next if ref($val);
                        next if $key =~ /^irods/;
                        $single_sample_meta->{$key} = $val;
                    }
                    $single_sample_meta->{sample} = $real_sample_name;
                    
                    my $dist_file = $self->output_file(sub_dir => $sub_dir, output_key => 'dist_file', basename => 'dist.dat', type => 'txt', metadata => $single_sample_meta);
                    my $plot_file = $self->output_file(sub_dir => $sub_dir, output_key => 'plot_file', basename => 'dist.py',  type => 'txt', metadata => $single_sample_meta);
                    
                    my @outfiles  = ($plot_file, $dist_file);
                    my $vcf_path  = $vcf->path->stringify;
                    my $dist_path = $dist_file->path;
                    
                    my @chroms = (1 .. 22, ("X", "Y", "MT"));
                    foreach my $chr (@chroms) {
                        my $png_file = $self->output_file(sub_dir => $sub_dir, output_key => 'png_files', basename => "dist.chr$chr.png", type => 'png', metadata => $single_sample_meta);
                        push(@outfiles, $png_file);
                        $self->relate_input_to_output($vcf_path, 'dist_plot', $png_file->path->stringify);
                    }
                    
                    my $this_cmd = "use VRPipe::Steps::polysomy; VRPipe::Steps::polysomy->run_and_check(vcf => q[$vcf_path], dist => q[$dist_path], bcftools => q[$bcftools_exe], bcftools_opts => q[$bcftools_opts], query => q[$query], sample => q[$real_sample_name], python => q[$python_exe]);";
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
            png_files => VRPipe::StepIODefinition->create(
                type            => 'png',
                min_files       => 0,
                max_files       => -1,
                description     => 'output plots by chr from polysomy',
                check_existence => 0,
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
    
    method run_and_check (ClassName|Object $self: Str|File :$vcf!, Str|File :$dist!, Str :$bcftools!, Str :$bcftools_opts!, Str :$query!, Str :$sample!, Str :$python!) {
        my $vcf_file  = VRPipe::File->get(path => $vcf);
        my $dist_file = VRPipe::File->get(path => $dist);
        my $dist_path = $dist_file->path->stringify;
        my $plot_path = $dist_path;
        $plot_path =~ s/\.dat$/.py/;
        
        my $cmd_line = "$bcftools polysomy " . $vcf_file->path . " $bcftools_opts -s $query -o " . $dist_file->dir;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $self->throw("File $dist doesn't exist..") unless -e $dist;
        
        # get the sample node in graph db
        my $vrtrack       = VRPipe::Schema->create('VRTrack');
        my $sample_source = $vrtrack->sample_source($sample);
        my $sample_props  = $vrtrack->sample_props_from_string($sample, $sample_source);
        my $sample_node   = $vrtrack->get('Sample', $sample_props);
        unless ($sample_node) {
            $self->throw("No sample node with name $sample_props->{name} was found in the graph db (for vcf $vcf, query $sample)");
        }
        
        # figure out the gender so we know what copy number to expect on chr X
        my %genders;
        foreach my $gender ($sample_node->related(outgoing => { type => 'gender', namespace => 'VRTrack', label => 'Gender' })) {
            $genders{ $gender->gender }++;
        }
        my @genders = keys %genders;
        my $gender = @genders == 1 ? $genders[0] : undef;
        
        # identify chrs with abnormal copy numbers, store them in metadata and
        # plot their BAF distribution
        $dist_file->update_stats_from_disc;
        my $fh = $dist_file->openr;
        my @ab_chrs;
        while (<$fh>) {
            next unless /^CN\s+(\S+)\s+(\S+)/;
            my $chr = $1;
            my $cn  = $2;
            
            my $expected = 2;
            if ($gender && $chr eq 'X') {
                if ($gender eq 'M') {
                    $expected = 1;
                }
            }
            
            if ($cn !~ /^$expected\.0+$/) {
                push(@ab_chrs, $chr);
            }
        }
        
        if (@ab_chrs) {
            my $chrs_str = $sample . ":" . join(",", @ab_chrs);
            $vcf_file->merge_metadata({ polysomy_chrs => $chrs_str });
            $sample_node->aberrant_chrs(\@ab_chrs);
            foreach my $chr (@ab_chrs) {
                my $cmd_line = "$python $plot_path -d $chr";
                system($cmd_line) && $self->throw("failed to run [$cmd_line]");
            }
        }
        
        $self->relate_input_to_output($vcf_file->path->stringify, 'polysomy_dist', $dist_path);
        $self->relate_input_to_output($vcf_file->path->stringify, 'polysomy_plot', $plot_path);
        
        return 1;
    }

}

1;
