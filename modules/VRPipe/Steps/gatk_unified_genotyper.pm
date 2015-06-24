
=head1 NAME

VRPipe::Steps::gatk_unified_genotyper - a step

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

#Example generic command for UnifiedGenotyper GATK v1.3
# java -jar GenomeAnalysisTK.jar \
#   -R resources/Homo_sapiens_assembly18.fasta \
#   -T UnifiedGenotyper \
#   -I sample1.bam [-I sample2.bam ...] \
#   -o snps.raw.vcf \
#   -stand_call_conf [50.0] \
#   -stand_emit_conf 10.0 \
#   -dcov [50] \

class VRPipe::Steps::gatk_unified_genotyper extends VRPipe::Steps::gatk {
    around options_definition {
        return {
            %{ $self->$orig }, # gatk options
            unified_genotyper_options => VRPipe::StepOption->create(description => 'Options for GATK UnifiedGenotyper, excluding -R,-I,-o'),
            minimum_records           => VRPipe::StepOption->create(description => 'Minimum number of records expected in output VCF. Not recommended if using genome chunking', optional => 1, default_value => 0)
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
            sites_file => VRPipe::StepIODefinition->create(type => 'vcf', min_files => 0, max_files => 1, description => 'Optional sites file for calling only at the given sites'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{bam_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta, caller => 'GATK_UnifiedGenotyper' };
            my $options = $self->handle_override_options($vcf_meta);
            $self->handle_standard_options($options);
            
            my $reference_fasta = $options->{reference_fasta};
            my $genotyper_opts  = $options->{unified_genotyper_options};
            my $minimum_records = $options->{minimum_records};
            
            if ($genotyper_opts =~ /$reference_fasta|-I|-o |UnifiedGenotyper/) {
                $self->throw("unified_genotyper_options should not include the reference, input or output options or UnifiedGenotyper task command");
            }
            
            my $summary_opts = $genotyper_opts;
            if ($self->inputs->{sites_file}) {
                $self->throw("unified_genotyper_options cannot contain the -alleles or --genotyping_mode (-gt_mode) options if a sites_file is an input to this step") if ($genotyper_opts =~ /-alleles/ || $genotyper_opts =~ /-gt_mode/ || $genotyper_opts =~ /--genotyping_mode/);
                my $sites_file = $self->inputs->{sites_file}[0];
                $summary_opts   .= ' --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles $sites_file';
                $genotyper_opts .= " --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles " . $sites_file->path;
                my $sites_meta = $sites_file->metadata;
                if (defined $$sites_meta{chrom} && defined $$sites_meta{from} && defined $$sites_meta{to}) {
                    if (defined $$vcf_meta{chrom} && defined $$vcf_meta{from} && defined $$vcf_meta{to}) {
                        unless ($$vcf_meta{chrom} eq $$sites_meta{chrom} && $$vcf_meta{from} == $$sites_meta{from} && $$vcf_meta{to} == $$sites_meta{to}) {
                            $self->throw("chrom/from/to metadata for output VCF and input sites file does not match");
                        }
                    }
                    else {
                        $$vcf_meta{chrom} = $$sites_meta{chrom};
                        $$vcf_meta{from}  = $$sites_meta{from};
                        $$vcf_meta{to}    = $$sites_meta{to};
                    }
                }
                if (defined $$sites_meta{seq_no}) {
                    ## need to add this so we can vcf-concat will work
                    if (defined $$vcf_meta{seq_no}) {
                        unless ($$vcf_meta{seq_no} eq $$sites_meta{seq_no}) {
                            $self->throw("seq_no metadata for output VCF and input sites file does not match");
                        }
                    }
                    else {
                        $$vcf_meta{seq_no} = $$sites_meta{seq_no};
                    }
                }
            }
            
            my $bams_list_path = $self->output_file(basename => 'bams.list', type => 'txt', temporary => 1)->path;
            my @input_ids = map { $_->id } @{ $self->inputs->{bam_files} };
            
            my $basename = 'gatk.vcf.gz';
            if (defined $$vcf_meta{chrom} && defined $$vcf_meta{from} && defined $$vcf_meta{to}) {
                my ($chrom, $from, $to) = ($$vcf_meta{chrom}, $$vcf_meta{from}, $$vcf_meta{to});
                $summary_opts   .= ' -L $region';
                $genotyper_opts .= " -L $chrom:$from-$to";
                $basename = "${chrom}_${from}-${to}.$basename";
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference_fasta -I $bams_list -o $vcf_file ' . $summary_opts
                )
            );
            
            my $vcf_file = $self->output_file(output_key => 'gatk_vcf_file', basename => $basename, type => 'vcf', metadata => $vcf_meta);
            my $vcf_path = $vcf_file->path;
            
            my $req      = $self->new_requirements(memory => 6000, time => 1);
            my $cmd      = $self->gatk_prefix($req->memory) . qq[ -T UnifiedGenotyper -R $reference_fasta -I $bams_list_path -o $vcf_path $genotyper_opts];
            my $this_cmd = "use VRPipe::Steps::gatk_unified_genotyper; VRPipe::Steps::gatk_unified_genotyper->genotype_and_check(q[$cmd], input_ids => [qw(@input_ids)], minimum_records => $minimum_records);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file] });
        };
    }
    
    method outputs_definition {
        return { gatk_vcf_file => VRPipe::StepIODefinition->create(type => 'vcf', max_files => 1, description => 'a single vcf file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run GATK UnifiedGenotyper for one or more BAMs, generating one compressed VCF per set of BAMs. Sites list can be provided";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method genotype_and_check (ClassName|Object $self: Str $cmd_line, ArrayRef[Int] :$input_ids!, Int :$minimum_records = 0) {
        my ($fofn_path, $out_path) = $cmd_line =~ /-I (\S+) -o (\S+)/;
        $fofn_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $fofn     = VRPipe::File->get(path => $fofn_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        my @input_files = map { VRPipe::File->get(id => $_) } @$input_ids;
        $fofn->create_fofn(\@input_files);
        
        $fofn->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        
        my $output_records = $out_file->num_records;
        if ($output_records < $minimum_records) {
            $out_file->unlink;
            $self->throw("Output VCF has $output_records data lines, fewer than required minimum $minimum_records");
        }
        else {
            return 1;
        }
    }
}

1;
