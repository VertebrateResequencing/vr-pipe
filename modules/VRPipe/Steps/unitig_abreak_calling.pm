
=head1 NAME

VRPipe::Steps::unitig_abreak_calling - a step

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

class VRPipe::Steps::unitig_abreak_calling with VRPipe::StepRole {
    method options_definition {
        return {
            htsbox_exe => VRPipe::StepOption->create(
                description   => 'path to htsbox executable',
                optional      => 1,
                default_value => 'htsbox'
            ),
            htsbox_abreak_options => VRPipe::StepOption->create(
                description   => 'options for htsbox abreak',
                optional      => 1,
                default_value => '-bcu'
            ),
            vcf_sort_exe => VRPipe::StepOption->create(
                description   => 'path to vcf-sort executable',
                optional      => 1,
                default_value => 'vcf-sort'
            ),
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to reference genome fasta'),
        };
    }
    
    method inputs_definition {
        return {
            name_sorted_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more name-sorted BAM files from which to call variants'
            ),
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options     = $self->options;
            my $htsbox      = $options->{htsbox_exe};
            my $abreak_opts = $options->{htsbox_abreak_options};
            my $vcf_sort    = $options->{vcf_sort_exe};
            
            my $reference_fasta = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $reference_fasta->is_absolute;
            
            if ($abreak_opts =~ /-f|$reference_fasta/) {
                $self->throw("htsbox_abreak_options should not include the reference or the -f option");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'htsbox',
                    version => VRPipe::StepCmdSummary->determine_version($htsbox, '^Version: (.+)$'),
                    summary => "htsbox abreak $abreak_opts -f \$reference_fasta \$name_sorted_bam_file | vcf-sort -c | bgzip -c > \$vcf_file && htsbox tabix -p vcf \$vcf_file"
                )
            );
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam (@{ $self->inputs->{name_sorted_bam_files} }) {
                my $prefix = $bam->basename;
                $prefix =~ s/\.bam$//;
                my $meta      = $bam->metadata;
                my $vcf_file  = $self->output_file(output_key => 'abreak_vcf_files', basename => "$prefix.vcf.gz", type => 'vcf', metadata => $meta);
                my $vcf_index = $self->output_file(output_key => 'abreak_vcf_index_files', basename => "$prefix.vcf.gz.tbi", type => 'bin', metadata => $meta);
                my $vcf_path  = $vcf_file->path;
                my $bam_path  = $bam->path;
                my $cmd       = qq[q[$htsbox abreak $abreak_opts -f $reference_fasta $bam_path | $vcf_sort -c | bgzip -c > $vcf_path && $htsbox tabix -p vcf $vcf_path]];
                my $this_cmd  = "use VRPipe::Steps::unitig_abreak_calling; VRPipe::Steps::unitig_abreak_calling->abreak_vcf_and_check($cmd);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$vcf_file, $vcf_index] });
            }
        };
    }
    
    method outputs_definition {
        return {
            abreak_vcf_files       => VRPipe::StepIODefinition->create(type => 'vcf', max_files => 1, description => 'a bgzipped .vcf.gz file for each set of input BAM files'),
            abreak_vcf_index_files => VRPipe::StepIODefinition->create(type => 'bin', max_files => 1, description => 'a tabix index file for each .vcf.gz file')
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Run htsbox abreak SV calling for one or more BAM files, generating one VCF file per BAM file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method abreak_vcf_and_check (ClassName|Object $self: Str $cmd_line!) {
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
