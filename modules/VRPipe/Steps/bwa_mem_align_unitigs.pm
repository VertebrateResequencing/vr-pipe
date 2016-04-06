
=head1 NAME

VRPipe::Steps::bwa_mem_align_unitigs - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::bwa_mem_align_unitigs extends VRPipe::Steps::bwa_mem_fastq {
    around options_definition {
        return {
            %{ $self->$orig },
            bwa_mem_options => VRPipe::StepOption->create(
                description   => 'options to bwa mem, excluding the input fastq, and reference',
                optional      => 1,
                default_value => '-x intractg -t 16'
            ),
            samtools_exe          => VRPipe::StepOption->create(description => 'path to your samtools executable',                                     optional => 1, default_value => 'samtools'),
            samtools_view_options => VRPipe::StepOption->create(description => 'samtools view option to control compression level of output bam file', optional => 1, default_value => '-1'),
            k8_exe                => VRPipe::StepOption->create(
                description   => 'path to your k8 executable',
                optional      => 1,
                default_value => 'k8'
            ),
            bwa_postalt => VRPipe::StepOption->create(
                description   => 'path to your bwa-postalt.js to be optionally used if mapping to a reference with ALTs',
                optional      => 1,
                default_value => ''
            ),
        };
    }
    
    around inputs_definition {
        return {
            %{ $self->$orig },
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                max_files   => -1,
                description => 'fastq files, which will be alnd independently',
                metadata    => {
                    lane           => 'lane name (a unique identifer for this sequencing run)',
                    library        => 'library name',
                    sample         => 'sample name',
                    center_name    => 'center name',
                    platform       => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study          => 'name of the study, put in the DS field of the RG header line',
                    analysis_group => 'project analysis group',
                    population     => 'sample population',
                    reads          => 'total number of reads (sequences)',
                    optional       => ['lane', 'library', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            my $ref     = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            my $k8      = $options->{k8_exe};
            my $postalt = $options->{bwa_postalt};
            
            my $bwa_exe  = $options->{bwa_exe};
            my $bwa_opts = $options->{bwa_mem_options};
            if ($bwa_opts =~ /$ref|mem/) {
                $self->throw("bwa_mem_options should not include the reference or the mem sub command");
            }
            my $cmd = $bwa_exe . ' mem ' . $bwa_opts;
            if ($self->inputs->{dict_file}) {
                $cmd .= ' -H ' . $self->inputs->{dict_file}->[0]->path;
            }
            my $samtools          = $options->{samtools_exe};
            my $samtools_compress = $options->{samtools_view_options};
            my $post_pipes        = qq[$samtools view $samtools_compress -];
            
            my $has_hla = 0;
            my $do_alt  = 0;
            if ($postalt && -f "$ref.alt") {
                $do_alt = 1;
                open(my $fh, "$ref.alt") || die("Could not open $ref.alt");
                while (<$fh>) {
                    $has_hla = 1 if /^HLA-[^\s\*]+\*\d+/;
                }
                close($fh) || die("Could not close $ref.alt");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bwa',
                    version => VRPipe::StepCmdSummary->determine_version($bwa_exe, '^Version: (.+)$'),
                    summary => 'bwa mem ' . $bwa_opts . ' $reference_fasta $fastq_file | samtools view ' . $samtools_compress . ' > $sam_file'
                )
            );
            
            my ($cpus) = $bwa_opts =~ m/-t\s*(\d+)/;
            my $req = $self->new_requirements(memory => 4900, time => 2, $cpus ? (cpus => $cpus) : ());
            
            my $idx = 0;
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                my $fq_meta = $fq->metadata;
                my $bam_meta = { reads => $fq_meta->{reads} };
                
                # add metadata and construct RG line
                my $rg_line;
                if (defined $fq_meta->{lane}) {
                    my $ln = $fq_meta->{lane};
                    $bam_meta->{lane} = $ln;
                    $rg_line = '@RG\tID:' . $ln;
                }
                else {
                    $rg_line = '@RG\tID:1';
                }
                if (defined $fq_meta->{library}) {
                    my $lb = $fq_meta->{library};
                    $bam_meta->{library} = $lb;
                    $rg_line .= '\tLB:' . $lb;
                }
                if (defined $fq_meta->{sample}) {
                    my $sm = $fq_meta->{sample};
                    $bam_meta->{sample} = $sm;
                    $rg_line .= '\tSM:' . $sm;
                }
                if (defined $fq_meta->{insert_size}) {
                    my $pi = $fq_meta->{insert_size};
                    $bam_meta->{insert_size} = $pi;
                    $rg_line .= '\tPI:' . $pi;
                }
                if (defined $fq_meta->{center_name}) {
                    my $cn = $fq_meta->{center_name};
                    $bam_meta->{center_name} = $cn;
                    $rg_line .= '\tCN:' . $cn;
                }
                if (defined $fq_meta->{platform}) {
                    my $pl = $fq_meta->{platform};
                    $bam_meta->{platform} = $pl;
                    $rg_line .= '\tPL:' . $pl;
                }
                if (defined $fq_meta->{study}) {
                    my $ds = $fq_meta->{study};
                    $bam_meta->{study} = $ds;
                    $rg_line .= '\tDS:' . $ds;
                }
                if (defined $fq_meta->{analysis_group}) {
                    $bam_meta->{analysis_group} = $fq_meta->{analysis_group};
                }
                if (defined $fq_meta->{population}) {
                    $bam_meta->{population} = $fq_meta->{population};
                }
                if (defined $fq_meta->{sample_id}) {
                    $bam_meta->{sample_id} = $fq_meta->{sample_id};
                }
                
                my $bam_file = $self->output_file(
                    output_key => 'bwa_mem_bam_files',
                    basename   => "$idx.unitig.bam",
                    type       => 'bam',
                    metadata   => $bam_meta
                );
                my @outfiles         = ($bam_file);
                my $these_post_pipes = $post_pipes;
                if ($do_alt) {
                    my $hla_pre = '';
                    if ($has_hla) {
                        my $fq_meta = $bam_meta;
                        delete $$fq_meta{reads} if exists $$fq_meta{reads};
                        delete $$fq_meta{bases} if exists $$fq_meta{bases};
                        foreach my $suffix (qw(HLA-A.fq HLA-B.fq HLA-C.fq HLA-DQA1.fq HLA-DQB1.fq HLA-DRB1.fq)) {
                            my $hla = $self->output_file(
                                sub_dir    => "$idx",
                                output_key => 'hla_fastq_files',
                                basename   => "$idx.unitig.hla.$suffix",
                                type       => 'fq',
                                metadata   => $fq_meta
                            );
                            push(@outfiles, $hla);
                        }
                        $hla_pre = ' -p ' . file($bam_file->dir, "$idx.unitig.hla");
                    }
                    $these_post_pipes = "$k8 $postalt$hla_pre $ref.alt | $these_post_pipes";
                }
                $idx++;
                
                my $this_cmd = "$cmd -R '$rg_line' $ref " . $fq->path . " | $these_post_pipes > " . $bam_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bwa_mem_align_unitigs', 'mem_and_check', [$this_cmd, $req, { output_files => \@outfiles }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            bwa_mem_bam_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => 'output bam file(s) mapped by bwa mem',
                metadata    => {
                    library        => 'library name',
                    sample         => 'sample name',
                    center_name    => 'center name',
                    platform       => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study          => 'name of the study, put in the DS field of the RG header line',
                    insert_size    => 'expected (mean) insert size if paired',
                    analysis_group => 'project analysis group',
                    population     => 'sample population',
                    reads          => 'total number of reads (sequences)',
                    optional       => ['library', 'insert_size', 'analysis_group', 'population', 'sample', 'center_name', 'platform', 'study']
                }
            ),
            hla_fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                min_files   => 0,
                max_files   => -1,
                description => 'optional output HLA fastq file(s) mapped by bwa mem',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Aligns the input unitig fastq with bwa mem to the reference; Produces an unsorted BAM file for each input fastq";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }

}

1;
