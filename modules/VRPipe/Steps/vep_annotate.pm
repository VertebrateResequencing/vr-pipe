
=head1 NAME

VRPipe::Steps::vep_annotate - a step

=head1 DESCRIPTION

Runs Ensembl's variant effect predictor to annotate input vcf files.

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

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

class VRPipe::Steps::vep_annotate with VRPipe::StepRole {
    method options_definition {
        return {
            vep_options       => VRPipe::StepOption->create(description => 'options to vep, excluding -i or -o',                                                                                                               optional      => 1),
            vep_exe           => VRPipe::StepOption->create(description => 'path to your vep script',                                                                                                                          optional      => 1, default_value => 'variant_effect_predictor.pl'),
            ensembl_api_paths => VRPipe::StepOption->create(description => 'paths to Ensembl Perl API installation directories separated by colon (Will to be prepended to PERL5LIB environment variable before running VEP)', optional      => 1),
            post_vep_options  => VRPipe::StepOption->create(description => 'after running VEP, option to pipe output vcf through another command (must contain either $output_vcf or $output_bcf string)',                     optional      => 1),
            bcftools_exe      => VRPipe::StepOption->create(description => 'path to bcftools executable (Will also be used to replace $bcftools placeholder in post_vep_options)',                                             default_value => 'bcftools'
            ),
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => '1 or more vcf or bcf files to annotate'),
        };
    }
    
    method body_sub {
        return sub {
            my $self              = shift;
            my $options           = $self->options;
            my $vep_exe           = $options->{vep_exe};
            my $vep_opts          = $options->{vep_options};
            my $bcftools_exe      = $options->{bcftools_exe};
            my $post_vep_cmds     = $options->{post_vep_options};
            my $ensembl_api_paths = $options->{ensembl_api_paths};
            if ($vep_opts =~ /-[i,o] /) {
                $self->throw("vep_options should not include the -i or -o option");
            }
            if (defined $post_vep_cmds && $post_vep_cmds !~ /\$output_vcf/ && $post_vep_cmds !~ /\$output_bcf/) {
                $self->throw("post_vep_options must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'variant_effect_predictor.pl',
                    version => VRPipe::StepCmdSummary->determine_version($vep_exe, '^version (.+)$'),
                    summary => "perl variant_effect_predictor.pl $vep_opts -i \$input -o \$output"
                )
            );
            my $suffix = $post_vep_cmds =~ /\$output_bcf/ ? "bcf" : "vcf.gz";
            foreach my $input_file (@{ $self->inputs->{vcf_files} }) {
                my $input_meta = $input_file->metadata;
                my $basename   = $input_file->basename;
                $basename =~ s/(\.vcf\.gz$|\.bcf$)//;
                #$self->output_file(basename => "$basename.${suffix}_temp.vcf", type => 'vcf', temporary => 1);
                my $output_file  = $self->output_file(output_key => 'vep_annot_vcf_files',       basename => "$basename.$suffix",     type => 'var', metadata => $input_meta);
                my $output_index = $self->output_file(output_key => 'vep_annot_vcf_index_files', basename => "$basename.$suffix.csi", type => 'idx', metadata => $input_meta);
                my $input_path   = $input_file->path;
                my $output_path  = $output_file->path;
                my $post_vep     = $post_vep_cmds;
                $post_vep =~ s/\$output_(b|v)cf/$output_path/;
                my $req = $self->new_requirements(memory => 5000, time => 1);
                my $this_cmd = "use VRPipe::Steps::vep_annotate; VRPipe::Steps::vep_annotate->vep_annotate_and_check(input => q[$input_path], chunk => [], output => q[$output_path], bcftools => q[$bcftools_exe], vep => q[$vep_exe], vep_opts => q[$vep_opts], api_paths => q[$ensembl_api_paths] , post_vep => q[$post_vep]);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$output_file, $output_index] });
            }
        };
    }
    
    method outputs_definition {
        return {
            vep_annot_vcf_files       => VRPipe::StepIODefinition->create(type => 'var', max_files => -1, description => 'output .bcf or .vcf files'),
            vep_annot_vcf_index_files => VRPipe::StepIODefinition->create(type => 'idx', max_files => -1, description => 'output CSI index for VCF/BCF files')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run Ensembl VEP to annotate one or more vcf or bcf files.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method vep_annotate_and_check (ClassName|Object $self: Str|File :$input!, ArrayRef[Str] :$chunk?, Str|File :$output!, Str :$bcftools!, Str :$vep!, Str :$vep_opts!, Str :$api_paths!, Str :$post_vep!) {
        my $input_file = VRPipe::File->get(path => $input);
        my $input_recs;
        my $region_opts = "";
        if (defined $chunk) {
            my $chrom = $$chunk[0];
            my $from  = $$chunk[1];
            my $to    = $$chunk[2];
            $region_opts = "-r $chrom:$from-$to | awk '\$0~/^#/ || \$2<$to'";
            $input_recs  = `$bcftools view -H $input $region_opts | wc -l`;
            chomp($input_recs);
        }
        else {
            $input_recs = $input_file->num_records;
        }
        
        #for empty chunks, copy over the header lines and exit
        unless ($input_recs) {
            my $opt = $output =~ /\.bcf$/ ? "-Ob" : "-Oz";
            my $cmd = "$bcftools view -h $input $opt -o $output && $bcftools index -c $output";
            system($cmd) && $self->throw("failed to run [$cmd]");
            return 1;
        }
        
        $input_file->disconnect;
        my $cmd_line = "export PERL5LIB=$api_paths; " if (defined $api_paths);
        #VEP v81 is buggy, doesn't compress with bgzip and -o STDOUT doesn't work with -quiet, therefore we output to a temporary vcf
        $cmd_line .= "$bcftools view $input $region_opts | $vep --quiet --vcf --no_progress --force_overwrite $vep_opts -o ${output}_temp.vcf";
        $cmd_line .= $post_vep ? " && cat ${output}_temp.vcf | $post_vep" : " && cat ${output}_temp.vcf | bgzip -c > $output";
        $cmd_line .= " && rm ${output}_temp.vcf && $bcftools index -c $output";
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output);
        $output_file->update_stats_from_disc;
        my $output_recs = `$bcftools view -H $output $region_opts | wc -l`;
        chomp($output_recs);
        
        unless ($output_recs == $input_recs) {
            $output_file->unlink;
            $self->throw("VEP output has $output_recs lines, less than input vcf records $input_recs");
        }
        else {
            return 1;
        }
    }
}

1;
