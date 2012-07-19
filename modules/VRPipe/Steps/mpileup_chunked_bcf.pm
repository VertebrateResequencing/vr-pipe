
=head1 NAME

VRPipe::Steps::mpileup_chunked_bcf - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::mpileup_chunked_bcf extends VRPipe::Steps::mpileup_bcf {
    around options_definition {
        return { %{ $self->$orig }, chunk_override_options => VRPipe::StepOption->create(description => 'File defining mpileup_bcf options to be overridden for specific chunks. This option is required, but may point to an empty file.'), };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type        => 'bam',
                                                               max_files   => -1,
                                                               description => '1 or more bam files to call variants',
                                                               metadata    => {
                                                                             chrom          => 'chrom',
                                                                             split_sequence => 'the chromosomal split sequence from the bam_split_by_sequence step',
                                                                             optional       => ['chrom', 'split_sequence'] }),
                 chunked_regions_file => VRPipe::StepIODefinition->create(type => 'txt', max_files => 1, description => 'File of chromosome region chunks to run concurrently') };
    }
    
    method body_sub {
        return sub {
            my $self             = shift;
            my $options          = $self->options;
            my $samtools         = $options->{samtools_exe};
            my $mpileup_opts     = $options->{samtools_mpileup_options};
            my $reference_fasta  = $options->{reference_fasta};
            my $override_options = $options->{chunk_override_options};
            my $override         = do $override_options;
            
            my %meta;
            my $bam_list;
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_path = $bam->path;
                $bam_list .= "$bam_path ";
                
                my $bam_meta = $bam->metadata;
                foreach my $key (qw(sample study project population continent analysis_group split_sequence chrom)) {
                    if (exists $$bam_meta{$key}) {
                        $meta{$key}->{ $$bam_meta{$key} } = 1;
                    }
                }
            }
            # Only keep unique metadata
            my %bcf_meta;
            foreach my $key (keys %meta) {
                my @vals = keys %{ $meta{$key} };
                next unless @vals == 1;
                $bcf_meta{$key} = $vals[0];
            }
            my $chrom_select = $bcf_meta{chrom} || $bcf_meta{split_sequence} || '';
            $chrom_select =~ s/^chrom//;
            
            my $seq_no     = 0;
            my $chunk_file = $self->inputs->{chunked_regions_file}[0];
            my $cfh        = $chunk_file->openr;
            while (<$cfh>) {
                my ($chr, $from, $to, undef) = split /\t/;
                next if ($chrom_select && !($chr =~ /^([cC]hr|[cC]hrom)?$chrom_select$/));
                my $region                = "${chr}_${from}-${to}";
                my $chunk_samtools        = $samtools;
                my $chunk_mpileup_options = $mpileup_opts;
                if (exists $$override{$region}{samtools_exe}) {
                    $chunk_samtools = $$override{$region}{samtools_exe};
                }
                if (exists $$override{$region}{samtools_mpileup_options}) {
                    $chunk_mpileup_options = $$override{$region}{samtools_mpileup_options};
                }
                
                $bcf_meta{seq_no} = $seq_no;
                $bcf_meta{region} = $region;
                ++$seq_no;
                my $bcf_file = $self->output_file(output_key => 'bcf_files', basename => "$region.bcf", type => 'bin', metadata => \%bcf_meta);
                my $bcf_path = $bcf_file->path;
                
                my $req = $self->new_requirements(memory => 500, time => 1);
                my $cmd = qq[$chunk_samtools mpileup $chunk_mpileup_options -r $chr:$from-$to -f $reference_fasta $bam_list > $bcf_path];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::mpileup_bcf', 'mpileup_bcf_and_check', [$cmd, $req, { output_files => [$bcf_file] }]);
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe     => 'samtools',
                                                                  version => VRPipe::StepCmdSummary->determine_version($samtools, '^Version: (.+)$'),
                                                                  summary => "samtools mpileup $mpileup_opts -r \$region -f \$reference_fasta \$bam_files > \$bcf_file"));
        };
    }
    
    method outputs_definition {
        return { bcf_files => VRPipe::StepIODefinition->create(type        => 'bin',
                                                               max_files   => -1,
                                                               description => 'a .bcf file for each chunked region',
                                                               metadata    => {
                                                                             chrom          => 'chrom',
                                                                             split_sequence => 'the chromosomal split sequence from the bam_split_by_sequence step',
                                                                             seq_no         => 'a sequence number assigned by the split for reassembly in correct order',
                                                                             optional       => ['chrom', 'split_sequence'] }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run samtools mpileup for one or more bams, generating one bcf file per set of bams";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
