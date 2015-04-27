
=head1 NAME

VRPipe::Steps::vcf_merge_different_samples - a step

=head1 DESCRIPTION

Run bcftools merge to merge a set of VCFs each with a different sample,
creating a single output VCF with all the input samples.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2014 Genome Research Limited.

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

class VRPipe::Steps::vcf_merge_different_samples extends VRPipe::Steps::bcftools {
    around options_definition {
        return {
            %{ $self->$orig },
            bcftools_options => VRPipe::StepOption->create(
                description   => 'Options for the bcftools merge command. Does not support the --print-header option; -O b or u is not supported.',
                optional      => 1,
                default_value => '-O z'
            )
        };
    }
    
    method inputs_definition {
        return {
            vcf_files => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'compressed vcf files',
                max_files   => -1
            )
        };
    }
    
    method _build_smaller_recommended_requirements_override {
        return 0;
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            my $bcftools_exe = $options->{bcftools_exe};
            my $bcfopts      = $options->{bcftools_options};
            
            if ($bcfopts =~ /-O\s*[bu]/) {
                # we're hard-coded for vcf output, so we really don't support
                # these bcf options
                $self->throw("-O b or u are not supported - this step can only output vcfs");
            }
            
            my $merged_basename = 'merged.vcf.gz';
            my $merged_meta     = $self->combined_metadata($self->inputs->{vcf_files});
            $self->_merge($bcftools_exe, $bcfopts, $self->inputs->{vcf_files}, $merged_basename, 'merged_vcf', 'vcf', $merged_meta, 'vcf_index');
        };
    }
    
    method outputs_definition {
        return {
            merged_vcf => VRPipe::StepIODefinition->create(
                type        => 'vcf',
                description => 'a merged vcf file',
                max_files   => 1
            ),
            vcf_index => VRPipe::StepIODefinition->create(
                type        => 'bin',
                description => 'index of the merged vcf file',
                max_files   => 1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Merges compressed VCFs using bcftools merge which contain different samples to produce a single VCF containing all the input samples";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method merge_vcf (ClassName|Object $self: Str $cmd_line) {
        my ($last_input_path, $output_path) = $cmd_line =~ /(\S+) > (\S+)/;
        
        my $last_input_file = VRPipe::File->get(path => $last_input_path);
        if ($last_input_file->type eq 'txt') {
            # it's a fofn, get the path of the first vcf in it
            my $fh   = $last_input_file->openr;
            my $line = <$fh>;
            $last_input_file->close;
            chomp($line);
            $last_input_file = VRPipe::File->get(path => $line);
        }
        my $input_lines = $last_input_file->lines;
        
        $last_input_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_lines = $output_file->lines;
        
        if ($output_lines < $input_lines) {
            $output_file->unlink;
            $self->throw("Output VCF has $output_lines, fewer than last input $input_lines");
        }
        else {
            return 1;
        }
    }
    
    method _merge (Str $bcftools_exe, Str $bcfopts, ArrayRef $input_vcfs, Str $output_basename, Str $output_key, Str $output_type, HashRef $output_metadata, Str $output_index_key) {
        # bcftools indexes and buffers all the input files in memory, amounting
        # to about 13.3MB memory usage each (including overhead), so we can get
        # a better memory estimate than VRPipe can guess.
        my $num_vcfs = scalar(@$input_vcfs);
        my $memory   = int($num_vcfs * 13.3);
        my $req      = $self->new_requirements(memory => $memory, time => 1);
        
        my $merged_vcf  = $self->output_file(output_key => $output_key, basename => $output_basename, type => $output_type, metadata => $output_metadata);
        my $output_path = $merged_vcf->path;
        my $index       = $self->output_file(output_key => $output_index_key, basename => $output_basename . '.csi', type => 'bin');
        
        my ($cmd_input, $summary_input);
        if ($num_vcfs == 1) {
            # merge doesn't work on 1 input file; just copy the file over.
            # (we don't symlink since it needs its own metadata, and we redo
            #  the metadata since copy adds source metadata to destination)
            my $source = $input_vcfs->[0];
            $source->copy($merged_vcf);
            $merged_vcf->add_metadata($output_metadata);
            $self->relate_input_to_output($source->path->stringify, 'copied', $output_path->stringify);
            $self->dispatch(["$bcftools_exe index $output_path", $req, { output_files => [$merged_vcf, $index] }]);
        }
        else {
            my ($cmd_input, $summary_input);
            my @input_paths = map { $_->path } @$input_vcfs;
            if ($num_vcfs < 27) {
                $cmd_input     = "@input_paths";
                $summary_input = '\@input_vcfs';
            }
            else {
                my $l_file = $self->output_file(temporary => 1, basename => $output_basename . '.l', type => 'txt');
                my $ofh = $l_file->openw;
                foreach my $path (@input_paths) {
                    print $ofh $path, "\n";
                }
                $l_file->close;
                $cmd_input     = '-l ' . $l_file->path;
                $summary_input = '-l $input_vcfs_fofn';
            }
            
            my $cmd = qq[$bcftools_exe merge $bcfopts $cmd_input > $output_path && $bcftools_exe index $output_path];
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => VRPipe::StepCmdSummary->determine_version($bcftools_exe, '^Version: (.+)$'),
                    summary => "bcftools merge $bcfopts $summary_input > \$output_path && bcftools index \$output_path"
                )
            );
            
            # relate input to output in graph db
            my @input_strings = map { $_->stringify } @input_paths;
            $self->relate_input_to_output(\@input_strings, 'merged', $output_path->stringify);
            
            $self->dispatch_wrapped_cmd(ref($self), 'merge_vcf', [$cmd, $req, { output_files => [$merged_vcf, $index] }]);
        }
    }
}

1;
