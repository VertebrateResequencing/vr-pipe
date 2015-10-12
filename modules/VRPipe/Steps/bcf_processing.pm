
=head1 NAME

VRPipe::Steps::bcf_processing - a step

=head1 DESCRIPTION

This is a generic vcf/bcf processing/filtering step. It takes bcf/vcf files  as
input, processes them using command(s) given in the options and outputs 
bcf/vcf files.

Example cmd_lines:  1) $bcftools view -Uf PASS -Ob -o $output_bcf $input_vcf 2)
bcftools norm -f ref.fa -Ou $input_bcf | bcftools annotate -a annots.tab.gz -h
annots.hdr -c CHROM,POS,REF,ALT,-,TAG -Ob > $output_bcf

Multiple pipes, semi-colons or redirections are allowed as long as one output 
vcf is produced for each input vcf. It is up to the user to provide a working 
command_line that produces correct outputs. The only requirement is that
strings  $input_vcf and $output_vcf are both specified in the command line.

Note the step should be used with caution as it does not carry out extra checks
 on the outputs as long as they are in bcf/vcf format.

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

class VRPipe::Steps::bcf_processing with VRPipe::StepRole {
    method options_definition {
        return {
            command_line => VRPipe::StepOption->create(
                description => 'Unix command(s) that will be used to produce the $output_vcf or $output_bcf from $input_vcf or $input_bcf (E.g. $bcftools view -Uf PASS -Ou $input_vcf | $bcftools annotate -x ID,FORMAT/DP -Ob -o $output_bcf)',
            ),
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to bcftools executable if using $bcftools placeholder in command_line',
                optional      => 1,
                default_value => 'bcftools'
            ),
            input_metadata_to_keep => VRPipe::StepOption->create(
                description => 'comma-separated list of metadata to be carried over to the processed files. By default copies all metadata to the output files; set to "-" to transfer no metadata',
                optional    => 1
            ),
            index_output => VRPipe::StepOption->create(
                description   => 'boolean; index the output VCF or BCF file when true',
                optional      => 1,
                default_value => 0
            ),
            check_records_vs_input => VRPipe::StepOption->create(
                description   => 'boolean; check the number of input records equals the number of output records when true',
                optional      => 1,
                default_value => 0
            ),
        };
    }
    
    method inputs_definition {
        return {
            bcf_files => VRPipe::StepIODefinition->create(
                type        => 'var',
                max_files   => -1,
                description => 'input VCF or BCF files to be processed',
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $cmd_line      = $options->{command_line};
            my $bcftools      = $options->{bcftools_exe};
            my $idx_output    = $options->{index_output};
            my $check_records = $options->{check_records_vs_input};
            my @metadata_keys = split(',', $options->{input_metadata_to_keep});
            
            if ($cmd_line !~ /\$input_vcf/ && $cmd_line !~ /\$input_bcf/) {
                $self->throw("cmd_line must contain one of the input strings \$input_vcf or \$input_bcf");
            }
            if ($cmd_line !~ /\$output_vcf/ && $cmd_line !~ /\$output_bcf/) {
                $self->throw("cmd_line must contain one of the output strings \$output_vcf or \$output_bcf");
            }
            if ($cmd_line =~ /\$bcftools/) {
                $self->set_cmd_summary(
                    VRPipe::StepCmdSummary->create(
                        exe     => 'bcftools',
                        version => VRPipe::StepCmdSummary->determine_version($bcftools, '^Version: (.+)$'),
                        summary => '$cmd_line'
                    )
                );
                $cmd_line =~ s/\$bcftools/$bcftools/g;
            }
            else {
                $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => '', version => '', summary => q[$cmd_line]));
            }
            
            my $sub_dir = 'a';
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $in_vcf (@{ $self->inputs->{bcf_files} }) {
                my $out_meta = {};
                my $in_meta  = $in_vcf->metadata;
                if (@metadata_keys) {
                    foreach my $key (@metadata_keys) {
                        next if $key eq '-';
                        $out_meta->{$key} = $in_meta->{$key};
                    }
                }
                else {
                    $out_meta = $in_meta;
                }
                $out_meta->{source_vcf} = $in_vcf->path;
                my $basename = $in_vcf->basename;
                $basename =~ s/\.[bv]cf(\.gz)?$//;
                my $suffix = $cmd_line =~ /\$output_vcf/ ? 'vcf.gz' : 'bcf';
                my $out_vcf     = $self->output_file(sub_dir => $sub_dir, output_key => 'processed_vcf_files', basename => "$basename.$suffix", type => $suffix, metadata => $out_meta);
                my @out_files   = ($out_vcf);
                my $input_path  = $in_vcf->path;
                my $output_path = $out_vcf->path;
                my $output_dir  = $out_vcf->dir;
                my $this_cmd    = $cmd_line;
                $this_cmd =~ s/\$input_[bv]cf/$input_path/g;
                $this_cmd =~ s/\$output_[bv]cf/$output_path/g;
                
                if ($idx_output) {
                    my $out_idx = $self->output_file(sub_dir => $sub_dir, output_key => 'processed_index_files', basename => "$basename.$suffix.csi", type => 'idx', metadata => $out_meta);
                    push @out_files, $out_idx;
                }
                $this_cmd = "cd $output_dir; " . $this_cmd;
                
                my $args = qq[q[$this_cmd], input_path => q[$input_path], output_path => q[$output_path]];
                $args .= qq[, check_records_vs_input => 1]  if $check_records;
                $args .= qq[, bcftools_exe => q[$bcftools]] if $idx_output;
                
                my $this_cmd_line = "use VRPipe::Steps::bcf_processing; VRPipe::Steps::bcf_processing->process_and_check($args);";
                $self->dispatch_vrpipecode($this_cmd_line, $req, { output_files => \@out_files });
                
                ++$sub_dir;
            }
        };
    }
    
    method outputs_definition {
        return {
            processed_vcf_files => VRPipe::StepIODefinition->create(
                type        => 'var',
                max_files   => -1,
                description => 'output VCF or BCF files produced by ad-hoc command',
                metadata    => { source_vcf => 'original BCF or VCF from which this file was generated' }
            ),
            processed_index_files => VRPipe::StepIODefinition->create(
                type        => 'idx',
                min_files   => 0,
                max_files   => -1,
                description => 'output CSI index for VCF/BCF files',
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Runs the input VCF or BCF through command line(s) given by the user";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method process_and_check (ClassName|Object $self: Str $cmd_line, Str|File :$input_path!, Str|File :$output_path!, Bool :$check_records_vs_input?, Str :$bcftools_exe?) {
        my $input_file  = VRPipe::File->get(path => $input_path);
        my $output_file = VRPipe::File->get(path => $output_path);
        
        $input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        if ($bcftools_exe) {
            system(qq[$bcftools_exe index -c $output_path]) && $self->throw("failed to index the output file [$bcftools_exe index -c $output_path]");
        }
        
        if ($check_records_vs_input) {
            $output_file->_filetype->check_records_vs_input($input_file, $cmd_line);
        }
    }

}

1;
