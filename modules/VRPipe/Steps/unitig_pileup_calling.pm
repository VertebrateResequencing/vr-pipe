
=head1 NAME

VRPipe::Steps::unitig_pileup_calling - a step

=head1 DESCRIPTION

...

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

class VRPipe::Steps::unitig_pileup_calling with VRPipe::StepRole {
    method options_definition {
        return {
            htsbox_exe => VRPipe::StepOption->create(
                description   => 'path to htsbox executable',
                optional      => 1,
                default_value => 'htsbox'
            ),
            htsbox_pileup_options => VRPipe::StepOption->create(
                description   => 'options for htsbox pileup',
                optional      => 1,
                default_value => '-Ccu'
            ),
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to reference genome fasta'),
        };
    }
    
    method inputs_definition {
        return {
            aln_files => VRPipe::StepIODefinition->create(
                type        => 'aln',
                max_files   => -1,
                description => '1 or more coordinate-sorted BAM files from which to call variants'
            ),
            sites_file => VRPipe::StepIODefinition->create(type => 'txt', min_files => 0, max_files => 1, description => 'Optional sites file for calling only at the given sites'),
        };
    }
    
    method body_sub {
        return sub {
            my $self     = shift;
            my $vcf_meta = $self->common_metadata($self->inputs->{aln_files});
            $vcf_meta = { %$vcf_meta, $self->element_meta };
            my $options = $self->handle_override_options($vcf_meta);
            
            my $htsbox      = $options->{htsbox_exe};
            my $pileup_opts = $options->{htsbox_pileup_options};
            
            my $reference_fasta = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $reference_fasta->is_absolute;
            
            if ($pileup_opts =~ /-f|-b|$reference_fasta/) {
                $self->throw("htsbox_pileup_options should not include the reference, the -f or the -b options");
            }
            
            if ($self->inputs->{sites_file}) {
                $self->throw("htsbox_pileup_options cannot contain the -b option if a sites_file is an input to this step") if ($pileup_opts =~ /-b/);
                my $sites_file = $self->inputs->{sites_file}->[0];
                $pileup_opts .= " -b " . $sites_file->path;
            }
            
            my @inputs = map { $_->path } @{ $self->inputs->{aln_files} };
            
            my $summary_opts = $pileup_opts;
            my $basename     = 'pileup.vcf.gz';
            if (defined $$vcf_meta{chrom} && defined $$vcf_meta{from} && defined $$vcf_meta{to}) {
                my ($chrom, $from, $to) = ($$vcf_meta{chrom}, $$vcf_meta{from}, $$vcf_meta{to});
                $summary_opts .= ' -r $region';
                $pileup_opts  .= " -r $chrom:$from-$to";
                $basename = "${chrom}_${from}-${to}.$basename";
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'htsbox',
                    version => VRPipe::StepCmdSummary->determine_version($htsbox, '^Version: (.+)$'),
                    summary => "htsbox pileup $summary_opts -f \$reference_fasta \$bam_file(s) | bgzip -c > \$vcf_file && htsbox tabix -p vcf \$vcf_file"
                )
            );
            
            my $vcf_file  = $self->output_file(output_key => 'unitig_pileup_vcf_files',       basename => $basename,       type => 'vcf', metadata => $vcf_meta);
            my $vcf_index = $self->output_file(output_key => 'unitig_pileup_vcf_index_files', basename => "$basename.tbi", type => 'bin', metadata => $vcf_meta);
            my $vcf_path  = $vcf_file->path;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            my $sample_id_fix = '';
            if (@inputs == 1 && $$vcf_meta{sample}) {
                my $input = $inputs[0];
                $input =~ s/\//\\\//g;
                my $sample_id = $$vcf_meta{sample};
                $sample_id_fix = qq[ | sed '/^#CHROM/s/$input\$/$sample_id/'];
            }
            my $cmd      = qq[q[$htsbox pileup $pileup_opts -f $reference_fasta @inputs$sample_id_fix | bgzip -c > $vcf_path && $htsbox tabix -p vcf $vcf_path]];
            my $this_cmd = "use VRPipe::Steps::unitig_pileup_calling; VRPipe::Steps::unitig_pileup_calling->pileup_vcf_and_check($cmd);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file, $vcf_index] });
        };
    }
    
    method outputs_definition {
        return {
            unitig_pileup_vcf_files       => VRPipe::StepIODefinition->create(type => 'vcf', max_files => 1, description => 'a bgzipped .vcf.gz file for each set of input BAM files'),
            unitig_pileup_vcf_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => 1, description => 'a tabix index file for each .vcf.gz file')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run htsbox pileup for one or more bams, generating one vcf file per set of bams";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method pileup_vcf_and_check (ClassName|Object $self: Str $cmd_line!) {
        my ($output_path) = $cmd_line =~ /> (\S+) &&/;
        $output_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $vcf = VRPipe::FileType->create('vcf', { file => $output_path });
        my $vcf_file = VRPipe::File->get(path => $output_path);
        
        $vcf_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $vcf_file->update_stats_from_disc(retries => 3);
        
        if ($vcf->check_type) {
            return 1;
        }
        else {
            $vcf_file->unlink;
            $self->throw("cmd [$cmd_line] failed because index file $output_path had the incorrect magic");
        }
    }
}

1;
