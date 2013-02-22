
=head1 NAME

VRPipe::Steps::gsnap - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

NJWalker <nw11@sanger.ac.uk>.

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

class VRPipe::Steps::gsnap with VRPipe::StepRole {
    use File::Basename;
    use Data::Dumper;
    
    method options_definition {
        return {
            gsnap_exe                            => VRPipe::StepOption->create(description => 'path to your gsnap executable',                                                                                                                              optional => 1, default_value => 'gsnap'),
            paired_end                           => VRPipe::StepOption->create(description => 'Set to 1 if input files are paired end. Default is for single end.',                                                                                         optional => 1, default_value => '0'),
            gsnap_db                             => VRPipe::StepOption->create(description => 'gsnap db that gsnap already knows about e.g. mm9, mm10, hg19, etc',                                                                                          optional => 1, default_value => 'mm9'),
            gsnap_genome_dir                     => VRPipe::StepOption->create(description => 'path to gsnap genome directory. If not set, will default to index known to the pipeline from a gmap_build step, or the default GMAP genome index directory', optional => 1),
            gsnap_input_gz                       => VRPipe::StepOption->create(description => 'gsnap input files are compressed with gunzip. Default true.',                                                                                                optional => 1, default_value => 1),
            gsnap_add_metadata_to_sam_read_group => VRPipe::StepOption->create(description => 'If they exist, values from metadata populate read-group id name lane and sample information in output sam',                                                  optional => 1, default_value => 1),
            
            readgroup_sm_from_metadata_key => VRPipe::StepOption->create(description => 'The SM of the readgroup will come from metadata associated with the fastq; this option chooses which metadata key to get the value from',                                                                                        optional => 1, default_value => 'sample'),
            gsnap_sam_use_0M               => VRPipe::StepOption->create(description => ' Insert 0M in CIGAR between adjacent insertions and deletions in output sam. Required by Picard, but can cause errors in other tools.',                                                                                          optional => 1, default_value => 1),
            gsnap_force_xs_dir             => VRPipe::StepOption->create(description => 'For RNA-Seq alignments, disallows XS:A:? when the sense direction is unclear, and replaces this value arbitrarily with XS:A:+. May be useful for some programs, such as Cufflinks, that cannot handle XS:A:?. See gsnap --help', optional => 1, default_value => 1),
            gsnap_format                   => VRPipe::StepOption->create(description => "Another format type, leave blank for gsnap default",                                                                                                                                                                             optional => 1, default_value => 'sam'),
            gsnap_other_options            => VRPipe::StepOption->create(description => "Other gsnap options",                                                                                                                                                                                                            optional => 1, default_value => '-t 12 -B 4 -N 1 --npaths=1 --filter-chastity=both --clip-overlap --fails-as-input --quality-protocol=sanger --format=sam')
        
        };
    }
    
    method inputs_definition {
        return {
            # sequence file - fastq for now
            fastq_files         => VRPipe::StepIODefinition->create(type => 'fq',  max_files => -1, description => '1 or more fastq files'),
            gmap_index_txt_file => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0,  max_files   => -1, description => '0 or 1 index file')
        };
    }
    
    method body_sub {
        return sub {
            my $self             = shift;
            my $options          = $self->options;
            my $gsnap_exe        = $options->{gsnap_exe};
            my $gsnap_db         = $options->{gsnap_db};
            my $paired           = $options->{paired_end};
            my $gsnap_genome_dir = "";
            if (defined($options->{gsnap_genome_dir})) {
                $gsnap_genome_dir = $options->{gsnap_genome_dir};
            }
            elsif (defined($self->inputs->{gmap_index_txt_file})) {
                my ($file) = @{ $self->inputs->{gmap_index_txt_file} };
                $gsnap_genome_dir = $file->dir->stringify;
            }
            else { $gsnap_genome_dir = undef; }
            my $gunzipped                            = $options->{gsnap_input_gz};
            my $gsnap_add_0M_to_cigar                = $options->{gsnap_sam_use_0M};
            my $gsnap_add_metadata_to_sam_read_group = $options->{gsnap_add_metadata_to_sam_read_group};
            my $sample_key                           = $options->{readgroup_sm_from_metadata_key};
            my $gsnap_force_xs_dir                   = $options->{gsnap_force_xs_dir};
            my $gsnap_format                         = $options->{gsnap_output};
            my $gsnap_other_options                  = $options->{gsnap_other_options};
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'gsnap', version => VRPipe::StepCmdSummary->determine_version($gsnap_exe . ' --version', 'GSNAP version  (.+) c'), summary => 'gsnap -d gsnap_db input_file'));
            my $req = $self->new_requirements(memory => 16000, time => 1); # more? 16GB RAM? Could be 8GB?
            my @input_file = @{ $self->inputs->{fastq_files} };
            my ($name) = fileparse($input_file[0]->basename, ('.fastq', '.fastq.gz', '.fq.gz', '.fq'));
            my ($cmd, $output_file_1, $output_file_2, $inputs);
            my $output_file_dir;
            my @outputfiles;
            
            # create command
            if ($paired) {
                $self->throw("Expecting two input files for paired end processing") unless @input_file == 2;
                $inputs = $input_file[0]->path . " " . $input_file[1]->path;
                
                $output_file_1 = $self->output_file(
                    output_key => 'gsnap_uniq_sam',
                    basename   => $name . ".concordant_uniq",
                    type       => 'txt',
                    metadata   => $input_file[0]->metadata
                );
                $output_file_dir = $output_file_1->dir->stringify;
            }
            
            if (!$paired) {
                # should only be one file
                $self->throw("Expecting one input file for single end processing") unless @input_file == 1;
                $inputs = $input_file[0]->path;
                
                $output_file_1 = $self->output_file(
                    output_key => 'gsnap_uniq_sam',
                    basename   => $name . ".unpaired_uniq",
                    type       => 'txt',
                    metadata   => $input_file[0]->metadata
                );
                $output_file_dir = $output_file_1->dir->stringify;
            }
            
            # deal with other options such as gunzip
            # construct command
            $cmd = "$gsnap_exe $inputs";
            $cmd .= " -D $gsnap_genome_dir "   if ($gsnap_genome_dir);
            $cmd .= " --gunzip "               if ($gunzipped);
            $cmd .= " --sam-use-0M "           if ($gsnap_add_0M_to_cigar);
            $cmd .= " --format=$gsnap_format " if ($gsnap_format);
            if ($gsnap_add_metadata_to_sam_read_group) {
                my $meta = $input_file[0]->metadata;
                my $library = $self->command_line_safe_string($meta->{library} || 'unknown_library');
                $cmd .= " --read-group-library $library";
                
                my $platform = $self->command_line_safe_string($meta->{platform} || 'unknown_platform');
                $cmd .= " --read-group-platform $platform";
                
                my $sample = $self->command_line_safe_string($meta->{$sample_key} || $meta->{sample} || 'unknown_sample');
                $cmd .= " --read-group-name $sample";
                
                my $study = $self->command_line_safe_string($meta->{study} || 'unknown_study');
                $cmd .= " --read-group-id $study";
            }
            $cmd .= " --force-xs-dir " if ($gsnap_force_xs_dir);
            $cmd .= " -d $gsnap_db " . $gsnap_other_options . " --split-output=$output_file_dir/$name";
            # -t 12 -B 4 -N 1 --npaths=1 --filter-chastity=both --clip-overlap --fails-as-input --quality-protocol=sanger --format=sam --split-output=$output_file_dir/$name";
            $self->dispatch([qq[$cmd], $req, { output_files => [$output_file_1] }]);
        };
    
    }
    
    method outputs_definition {
        return { gsnap_uniq_sam => VRPipe::StepIODefinition->create(type => 'txt', description => 'gsnap mapped sequences files in sam format'), };
    }
    
    method description {
        return "Step for GSNAP mapper";
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method command_line_safe_string (Str $str) {
        # add single quotes around the string if it contains spaces
        if ($str =~ /\s/) {
            $str = qq['$str'];
        }
        return $str;
    }
}

