
=head1 NAME

VRPipe::Steps::gatk_gvcf_linear_index - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>

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

class VRPipe::Steps::gatk_gvcf_linear_index extends VRPipe::Steps::gatk_v2 {
    around options_definition {
        return {
            %{ $self->$orig },
            gatk_cat_variants_options => VRPipe::StepOption->create(description => 'command line options for GATK CatVariants; excludes --variant, -R, and -out options which are set by another StepOption', optional => 1, default_value => '-assumeSorted')
        };
    }
    
    method inputs_definition {
        return {
            gvcf_files       => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => '1 or more gvcf files'),
            gvcf_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'index files for the input gvcf files'),
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            my $cat_opts = $options->{gatk_cat_variants_options};
            
            my $reference_fasta = $options->{reference_fasta};
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => qq[java \$jvm_args -jar -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R \$reference_fasta $cat_opts --variant \$gvcf.vcf.gz -out \$gvcf.g.vcf],
                )
            );
            
            my $req = $self->new_requirements(memory => 6000, time => 1);
            my $idx = 1;
            foreach my $gvcf (@{ $self->inputs->{gvcf_files} }) {
                my $gvcf_path = $gvcf->path;
                my $basename  = "$idx." . $gvcf->basename;
                my $meta      = { %{ $gvcf->metadata }, $self->element_meta };
                $basename =~ s/\.vcf\.gz$/.g.vcf/;
                my $uncompressed_gvcf_file  = $self->output_file(output_key => 'uncompressed_gvcf_files',       basename => $basename,          type => 'vcf', metadata => $meta);
                my $uncompressed_gvcf_index = $self->output_file(output_key => 'uncompressed_gvcf_index_files', basename => $basename . ".idx", type => 'bin', metadata => $meta);
                my $jvm_args                = $self->jvm_args($req->memory);
                my $cmd                     = $self->java_exe . qq[ $jvm_args -cp ] . $self->jar . qq[ org.broadinstitute.gatk.tools.CatVariants -R $reference_fasta $cat_opts --variant $gvcf_path -out ] . $uncompressed_gvcf_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::gatk_gvcf_linear_index', 'uncompress_and_check', [$cmd, $req, { output_files => [$uncompressed_gvcf_file, $uncompressed_gvcf_index] }]);
                $idx++;
            }
        };
    }
    
    method outputs_definition {
        return {
            uncompressed_gvcf_files       => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'uncompressed gvcf files'),
            uncompressed_gvcf_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => -1, description => 'linear index of uncompressed gvcf files'),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Decompress gVCF files and produce the linear index preferred by GATK for CombineGVCFs and GenotypeGVCFs.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method uncompress_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($input_path, $out_path) = $cmd_line =~ /--variant (\S+) -out (\S+)/;
        $input_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path   || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $input_file = VRPipe::File->get(path => $input_path);
        my $out_file   = VRPipe::File->get(path => $out_path);
        my $out_index  = VRPipe::File->get(path => $out_path . '.idx');
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc;
        $out_index->update_stats_from_disc;
        
        my $input_recs  = $input_file->num_records;
        my $output_recs = $out_file->num_records;
        
        unless ($output_recs == $input_recs) {
            $out_file->unlink;
            $out_index->unlink;
            $self->throw("Output gVCF has different number of data lines from input gVCF (input $input_recs, output $output_recs)");
        }
        else {
            return 1;
        }
    }
}

1;
