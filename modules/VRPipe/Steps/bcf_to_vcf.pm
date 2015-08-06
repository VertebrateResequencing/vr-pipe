
=head1 NAME

VRPipe::Steps::bcf_to_vcf - a step

=head1 DESCRIPTION

Generates a compressed VCF from a BCF file using bcftools call, optionally
restricting output on samples or sites

=head1 AUTHOR

Chris Joyce    <cj5@sanger.ac.uk>. Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2014 Genome Research Limited.

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

class VRPipe::Steps::bcf_to_vcf extends VRPipe::Steps::bcftools {
    around options_definition {
        return {
            %{ $self->$orig },
            bcftools_call_options => VRPipe::StepOption->create(
                description => 'bcftools calling options; v0 defaults to "-p 0.99 -vcgN"; v1 defaults to "-m"',
                optional    => 1
            ),
            sample_sex_file => VRPipe::StepOption->create(
                description => 'Tab- or space- delimited file listing each sample name and sex (M or F). If not provided, will call on all samples in the bcf header. If provided, calls will be made on the intersection of the samples in the file and samples in the bcf header.',
                optional    => 1
            ),
            assumed_sex => VRPipe::StepOption->create(
                description   => 'If M or F is not present for a sample in the sample sex file, then this sex is assumed',
                optional      => 1,
                default_value => 'F'
            ),
            minimum_records => VRPipe::StepOption->create(
                description   => 'Minimum number of records expected in output VCF. Not recommended if using genome chunking',
                optional      => 1,
                default_value => 0
            ),
            post_calling_vcftools => VRPipe::StepOption->create(
                description => 'After calling with bcftools call, option to pipe output vcf through a bcftools or other vcf streaming command, e.g. "$bcftools norm -f hs37d5.fa", where $bcftools will be replaced by the bcftools_exe option',
                optional    => 1
            ),
            vcf_sample_from_metadata => VRPipe::StepOption->create(
                description => 'if the sample id in the resulting vcf header matches metadata with key x, but you want it to match the value from key y, provide x:y; separate multiple y keys with + symbols - values will be joined with underscores. This only works with single-sample vcfs',
                optional    => 1
            )
        };
    }
    
    method inputs_definition {
        return {
            bcf_files  => VRPipe::StepIODefinition->create(type => 'bcf', max_files => -1, description => '1 or more bcf files to convert to compressed vcf'),
            sites_file => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0,  max_files   => 1, description => 'Optional sites file for calling only at the given sites'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{bcf_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($vcf_meta);
            
            my $bcftools  = $options->{bcftools_exe};
            my $call_opts = $options->{bcftools_call_options};
            $call_opts ||= $self->_bcftools_calling_defaults;
            my $calling_command = $self->_bcftools_calling_command;
            my $samples_option  = $self->_bcftools_samples_option;
            my $assumed_sex     = $options->{assumed_sex};
            my $minimum_records = $options->{minimum_records};
            my $sfm             = $options->{vcf_sample_from_metadata};
            my $post_filter     = $options->{post_calling_vcftools};
            $post_filter =~ s/\$bcftools/$bcftools/g;
            
            my $sample_sex_file;
            if ($options->{sample_sex_file}) {
                $sample_sex_file = file($options->{sample_sex_file});
                $self->throw("sample_sex_file must be an absolute path") unless $sample_sex_file->is_absolute;
            }
            if ($self->inputs->{sites_file}) {
                $self->throw("bcftools_call_options cannot contain the -l/-R/-T option if a sites_file is an input to this step") if ($call_opts =~ /-[lR]/); # can't throw on -T since that was a option with a different meaning in v0
                $call_opts .= $self->_bcftools_site_files_option($self->inputs->{sites_file}[0]->path);
            }
            
            my $output = $self->_bcftools_compressed_vcf_output($post_filter);
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bcf (@{ $self->inputs->{bcf_files} }) {
                my $bcf_meta = $bcf->metadata;
                my $bcf_path = $bcf->path;
                my $basename = $bcf->basename;
                $basename =~ s/\.bcf$//;
                $bcf_meta->{caller} = 'samtools_mpileup_bcftools';
                
                my $female_ploidy = defined $bcf_meta->{female_ploidy} ? delete $bcf_meta->{female_ploidy} : 2;
                my $male_ploidy   = defined $bcf_meta->{male_ploidy}   ? delete $bcf_meta->{male_ploidy}   : 2;
                
                my $temp_samples_path = $self->output_file(basename => $basename . '.samples', type => 'txt', temporary => 1)->path;
                
                my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => $basename . '.vcf.gz', type => 'vcf', metadata => $bcf_meta);
                my $vcf_path = $vcf_file->path;
                my $cmd_line = qq[$bcftools $calling_command $call_opts $samples_option $temp_samples_path $bcf_path $output > $vcf_path];
                
                my $bcf_id = $bcf->id;
                my $args   = qq['$cmd_line', '$temp_samples_path', source_file_ids => ['$bcf_id'], female_ploidy => '$female_ploidy', male_ploidy => '$male_ploidy', assumed_sex => '$assumed_sex'];
                $args .= qq[, sample_sex_file => '$sample_sex_file'] if $sample_sex_file;
                $args .= qq[, vcf_sample_from_metadata => '$sfm']    if $sfm;
                my $cmd = "use VRPipe::Steps::bcf_to_vcf; VRPipe::Steps::bcf_to_vcf->bcftools_call_with_sample_file($args, minimum_records => $minimum_records);";
                $self->dispatch_vrpipecode($cmd, $req, { output_files => [$vcf_file] });
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'bcftools',
                    version => $self->bcftools_version_string,
                    summary => "bcftools $calling_command $call_opts $samples_option \$samples_file \$bcf_file $output > \$vcf_file"
                )
            );
        };
    }
    
    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->create(type => 'vcf', max_files => -1, description => 'a .vcf.gz file for each input bcf file') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run bcftools call option to generate one compressed vcf file per input bcf. Will take care of ploidy if sample_sex_file is provided and bcf files contain male_ploidy and female_ploidy metadata.";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method bcftools_call_with_sample_file (ClassName|Object $self: Str $cmd_line!, Str|File $sample_ploidy_path, ArrayRef[Int] :$source_file_ids!, Str|File :$sample_sex_file?, Str :$vcf_sample_from_metadata?, Int :$female_ploidy!, Int :$male_ploidy!, Str :$assumed_sex = 'F', Int :$minimum_records = 0) {
        my @input_files = map { VRPipe::File->get(id => $_) } @$source_file_ids;
        
        # find out the samples contained in the input files
        my %source_samples;
        my $bam_input = 0;
        foreach my $file (@input_files) {
            if ($file->type eq 'bcf') {
                my $bcf_ft = VRPipe::FileType->create('bcf', { file => $file->path });
                map { $source_samples{$_} = 1 } @{ $bcf_ft->samples };
            }
            elsif ($file->type eq 'bam') {
                my $bp = VRPipe::Parser->create('bam', { file => $file });
                map { $source_samples{$_} = 1 } $bp->samples;
                $bp->close;
                $bam_input = 1;
            }
            else {
                $self->throw("Source file not of type bcf or bam, " . $file->path);
            }
        }
        $self->throw("No samples found") unless (keys %source_samples);
        
        # if the inputs to the step are bam files, this is the mpileup step
        # and we need to create a fofn of the input bam files
        if ($bam_input) {
            my ($bam_fofn_path) = $cmd_line =~ /-b (\S+) \|/;
            $self->throw("No bam fofn path found in command line $cmd_line") unless $bam_fofn_path;
            my $bam_fofn = VRPipe::File->get(path => $bam_fofn_path);
            $bam_fofn->create_fofn(\@input_files);
        }
        
        # find out the sex of each of the input files - either from
        # the sample_sex_file or assumed_sex
        # if sample_sex_file is a subset of all the samples, then
        # we'll be calling on that subset
        my %samples;
        if ($sample_sex_file) {
            my $sex_file = VRPipe::File->create(path => $sample_sex_file);
            my $fh = $sex_file->openr;
            while (<$fh>) {
                chomp;
                my ($sample, $sex) = split;
                next unless (exists $source_samples{$sample});
                $sex ||= $assumed_sex;
                $samples{$sample} = $sex;
            }
            $fh->close;
        }
        else {
            foreach my $sample (keys %source_samples) {
                $samples{$sample} = $assumed_sex;
            }
        }
        
        # write the samples and ploidy to a file which will be
        # used by bcftools to call
        my $ploidy_file      = VRPipe::File->get(path => $sample_ploidy_path);
        my $pfh              = $ploidy_file->openw;
        my $expected_samples = 0;
        while (my ($sample, $sex) = each %samples) {
            my $ploidy = $sex eq 'F' ? $female_ploidy : $male_ploidy;
            next unless $ploidy; ## ploidy may be zero. eg chromY for female samples
            print $pfh "$sample\t$ploidy\n";
            $expected_samples++;
        }
        $pfh->close;
        $ploidy_file->update_stats_from_disc(retries => 3);
        
        if ($sample_sex_file && !$expected_samples) {
            $self->throw("no samples in $sample_sex_file overlapped with samples in the input files");
        }
        
        # check that the number of samples is as expected
        my $actual_samples = $ploidy_file->lines;
        unless ($actual_samples == $expected_samples) {
            $ploidy_file->unlink;
            $self->throw("write_sample_ploidy_file failed because $actual_samples samples were generated in the output ploidy files, yet we expected $expected_samples samples");
        }
        
        # run the command line
        my ($output_path) = $cmd_line =~ /> (\S+)$/;
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        $output_file->update_stats_from_disc(retries => 3);
        
        my $ft = VRPipe::FileType->create('vcf', { file => $output_path });
        unless ($ft->num_header_lines > 0) {
            $output_file->unlink;
            $self->throw("Output VCF [$output_path] has no header lines");
        }
        
        # reheader the output vcf if we need to alter sample identifiers
        #*** this is a little gross since we're overwriting our existing output
        # file and haven't specified the needed temp files, and haven't asked
        # for bcftools and tabix exes, nor published that we run these commands
        # anywhere - this should probably happen in a dedicated step, but we
        # need a quick hack now
        if ($vcf_sample_from_metadata) {
            my $meta = $output_file->metadata;
            my ($src_key, $dst_key) = split(':', $vcf_sample_from_metadata);
            my @dst_keys = split(/\+/, $dst_key);
            
            if ($src_key && @dst_keys && defined $meta->{$src_key} && defined $meta->{ $dst_keys[0] }) {
                my $dir          = $output_file->dir;
                my $header_file  = VRPipe::File->create(path => file($dir, $output_file->basename . '.temp_header'));
                my $reheader_vcf = VRPipe::File->create(path => file($dir, $output_file->basename . '.temp_reheader.vcf.gz'), metadata => $output_file->metadata);
                
                my $cmd = "bcftools view -h $output_path";
                $header_file->disconnect;
                open(my $bvfh, "$cmd |") || $self->throw("Couldn't open pipe from [$cmd]\n");
                my $ofh = $header_file->openw;
                while (<$bvfh>) {
                    if (/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t([^\t]*\S)$/) {
                        my $sample = $1;
                        if ($meta->{$src_key} eq $sample) {
                            my $new_sample = join('_', map { $meta->{$_} || 'undef' } @dst_keys);
                            $_ =~ s/\t$sample/\t$new_sample/;
                        }
                    }
                    print $ofh $_;
                }
                $header_file->close || $self->throw("Couldn't close pipe from [$cmd]\n");
                
                $cmd = "tabix -r " . $header_file->path . " $output_path > " . $reheader_vcf->path;
                system($cmd) && $self->throw("Failed to run [$cmd]\n");
                
                $header_file->rm;
                $reheader_vcf->update_stats_from_disc;
                $reheader_vcf->mv($output_file);
                $output_file->update_stats_from_disc;
            }
        }
        
        my $output_records = $output_file->num_records;
        if ($output_records < $minimum_records) {
            $output_file->unlink;
            $self->throw("Output VCF has $output_records data lines, fewer than required minimum $minimum_records");
        }
        else {
            return 1;
        }
    }
}

1;
