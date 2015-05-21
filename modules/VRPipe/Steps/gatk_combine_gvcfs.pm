
=head1 NAME

VRPipe::Steps::gatk_combine_gvcfs - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk> and Yasin Memari <ym3@sanger.ac.uk>.

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

# java -Xmx2g -jar GenomeAnalysisTK.jar \
#   -R ref.fasta \
#   -T CombineGVCFs \
#   --variant gvcf1.vcf \
#   --variant gvcf2.vcf \
#   -o mergeGvcf.vcf

class VRPipe::Steps::gatk_combine_gvcfs extends VRPipe::Steps::gatk_v2 {
    around options_definition {
        return {
            %{ $self->$orig }, # gatk options
            combine_gvcfs_options    => VRPipe::StepOption->create(description => 'Options for GATK CombineGVCFs, excluding -R,-V,-o', optional      => 1),
            maximum_gvcfs_to_combine => VRPipe::StepOption->create(description => 'The maximum gVCFs to combine at once',              default_value => 200),
        };
    }
    
    method inputs_definition {
        return {
            gvcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => '1 or more gvcf files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{gvcf_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($vcf_meta);
            $self->handle_standard_options($options);
            
            my $reference_fasta          = $options->{reference_fasta};
            my $combine_gvcfs_opts       = $options->{combine_gvcfs_options} ? $options->{combine_gvcfs_options} : "";
            my $maximum_gvcfs_to_combine = $options->{maximum_gvcfs_to_combine};
            
            if ($combine_gvcfs_opts =~ /$reference_fasta|-V |--variant|-o | --output|CombineGVCFs/) {
                $self->throw("combine_gvcfs_options should not include the reference, input or output options or CombineGVCFs task command");
            }
            
            my $summary_opts = $combine_gvcfs_opts;
            my $basename     = 'gatk_mergeGvcf.vcf.gz';
            if (defined $$vcf_meta{chrom} && defined $$vcf_meta{from} && defined $$vcf_meta{to}) {
                my ($chrom, $from, $to) = ($$vcf_meta{chrom}, $$vcf_meta{from}, $$vcf_meta{to});
                $summary_opts       .= ' -L $region';
                $combine_gvcfs_opts .= " -L $chrom:$from-$to";
                $basename = "${chrom}_${from}-${to}.$basename";
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T CombineGVCFs -R $reference_fasta --variant $gvcf1 --variant $gvcf2 [...] -o $out_gvcf ' . $summary_opts
                )
            );
            
            my @input_files = @{ $self->inputs->{gvcf_files} };
            
            my $count = 0;
            while (@input_files) {
                my @files_list = splice @input_files, 0, $maximum_gvcfs_to_combine;
                my @file_inputs = map { "--variant " . $_->path } @files_list;
                my $vcf_file = $self->output_file(output_key => 'combine_gvcf_files', basename => "batch_$count." . $basename, type => 'vcf', metadata => $vcf_meta);
                my $vcf_path = $vcf_file->path;
                $self->output_file(output_key => 'vcf_index', basename => "batch_$count." . $basename . ".tbi", type => 'bin', metadata => $vcf_meta);
                
                my $req      = $self->new_requirements(memory => 6000, time => 1);
                my $jvm_args = $self->jvm_args($req->memory);
                my $cmd      = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T CombineGVCFs -R $reference_fasta @file_inputs -o $vcf_path $combine_gvcfs_opts];
                my $this_cmd = "use VRPipe::Steps::gatk_combine_gvcfs; VRPipe::Steps::gatk_combine_gvcfs->genotype_and_check(q[$cmd]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file] });
                $count++;
            }
        };
    }
    
    method outputs_definition {
        return {
            combine_gvcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'combined gvcf files'),
            vcf_index          => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'index of combined gvcf files'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run GATK CombineGVCFs on equal sized batches of input gVCF files produced by the Haplotype Caller, generating one joint VCF file for each batch";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method genotype_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $out_path) = $cmd_line =~ /--variant (\S+) .*? -o (\S+)/;
        $input_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path   || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $input_file = VRPipe::File->get(path => $input_path);
        my $out_file   = VRPipe::File->get(path => $out_path);
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        return 1;
    }
}

1;
