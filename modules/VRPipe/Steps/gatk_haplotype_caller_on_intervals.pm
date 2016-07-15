
=head1 NAME

VRPipe::Steps::gatk_haplotype_caller_on_intervals - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2016 Genome Research Limited.

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

class VRPipe::Steps::gatk_haplotype_caller_on_intervals extends VRPipe::Steps::gatk_haplotype_caller {
    around options_definition {
        return {
            %{ $self->$orig },
            list_of_intervals_files => VRPipe::StepOption->create(description => 'path to a file containing a list of interval list files. HaplotypeCaller will be run for each interval file in this list'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{bam_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta, caller => 'GATK_HaplotypeCaller' };
            my $options = $self->handle_override_options($vcf_meta);
            $self->handle_standard_options($options);
            
            my $reference_fasta = $options->{reference_fasta};
            my $haplotyper_opts = $options->{haplotype_caller_options};
            my $contam_key      = $options->{contamination_metadata_key};
            my $minimum_records = $options->{minimum_records};
            my $tabix           = $options->{tabix_exe};
            my $avx             = $options->{avx_lsf_requirement_string};
            my $intervals       = $options->{list_of_intervals_files};
            
            if ($haplotyper_opts =~ /$reference_fasta|-I |--input_file|-o | --output|HaplotypeCaller/) {
                $self->throw("haplotype_caller_options should not include the reference, input or output options or HaplotypeCaller task command");
            }
            
            my $summary_opts = qq[java \$jvm_args -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R \$reference_fasta -I \$bams_list -o \$vcf_file $haplotyper_opts -L \$region];
            if ($self->inputs->{bam_recalibration_files}) {
                my $bqsr_file = $self->inputs->{bam_recalibration_files}[0];
                $summary_opts .= ' -BQSR $recal_file';
                $haplotyper_opts .= ' --BQSR ' . $bqsr_file->path;
            }
            if ($contam_key && exists $$vcf_meta{$contam_key}) {
                $summary_opts    .= qq[ -contamination $contam_key];
                $haplotyper_opts .= qq[ -contamination $$vcf_meta{$contam_key}];
            }
            
            my $file_list_id;
            if (@{ $self->inputs->{bam_files} } > 1) {
                $file_list_id = VRPipe::FileList->create(files => $self->inputs->{bam_files})->id;
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => $summary_opts
                )
            );
            
            my ($cpus) = $haplotyper_opts =~ m/-nct\s*(\d+)/;
            unless ($cpus) {
                ($cpus) = $haplotyper_opts =~ m/--num_cpu_threads_per_data_thread\s*(\d+)/;
            }
            my $req = $self->new_requirements(memory => 8000, time => 1, $cpus ? (cpus => $cpus) : (), $avx ? (custom => { lsf => $avx }) : ());
            my $basename = 'gatk_haplotype.vcf.gz';
            
            my $seq_no = 0;
            open(my $fh, "<$intervals") || $self->throw("Could not open list_of_intervals_files $intervals");
            while (my $interval = <$fh>) {
                chomp $interval;
                next unless $interval;
                my $chunk_opts     = "$haplotyper_opts -L $interval";
                my $chunk_basename = "$seq_no.$basename";
                my $chunk_meta     = { %$vcf_meta, seq_no => $seq_no++ };
                
                my $bams_list_path;
                if ($file_list_id) {
                    $bams_list_path = $self->output_file(basename => "$chunk_basename.bams.list", type => 'txt', temporary => 1)->path;
                }
                else {
                    $bams_list_path = $self->inputs->{bam_files}->[0]->path;
                }
                
                my $vcf_file       = $self->output_file(output_key => 'gatk_hc_vcf_files',       basename => $chunk_basename,       type => 'vcf', metadata => $chunk_meta);
                my $vcf_index_file = $self->output_file(output_key => 'gatk_hc_vcf_index_files', basename => "$chunk_basename.tbi", type => 'tbi', metadata => $chunk_meta);
                my $vcf_path       = $vcf_file->path;
                
                $chunk_opts .= qq[ && $tabix -f -p vcf $vcf_path] if ($chunk_opts =~ m/--disable_auto_index_creation_and_locking_when_reading_rods/);
                my $cmd = 'q[' . $self->gatk_prefix($req->memory) . qq[ -T HaplotypeCaller -R $reference_fasta -I $bams_list_path -o $vcf_path $chunk_opts] . ']';
                $cmd .= qq[, input_file_list => $file_list_id] if $file_list_id;
                my $this_cmd = "use VRPipe::Steps::gatk_haplotype_caller; VRPipe::Steps::gatk_haplotype_caller->genotype_and_check($cmd, minimum_records => $minimum_records);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file, $vcf_index_file] });
            }
            close($fh) || $self->throw("Could not close list_of_intervals_files $intervals");
        };
    }
    
    method outputs_definition {
        return {
            gatk_hc_vcf_files       => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a single vcf file per interval file in the given list'),
            gatk_hc_vcf_index_files => VRPipe::StepIODefinition->create(type => 'tbi', max_files => -1, description => 'a single vcf index file per vcf'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run GATK HaplotypeCaller for one or more BAMs, generating one compressed VCF per set of BAMs for each of the given intervals files";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
