=head1 NAME

VRPipe::Steps::mpileup_vcf - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::mpileup_vcf extends VRPipe::Steps::mpileup_bcf {
    around options_definition {
	return { %{$self->$orig},
		 bcftools_exe => VRPipe::StepOption->get(description => 'path to bcftools executable',
							 optional => 1,
							 default_value => 'bcftools'),
		 bcftools_view_options => VRPipe::StepOption->get(description => 'bcftools view options',
								  optional => 1,
								  default_value => '-gcv'),
                 mimimum_calls => VRPipe::StepOption->get(description => 'minumum expected number of variant calls',
							  optional => 1,
							  default_value => 0),};
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $samtools = $options->{samtools_exe};
            my $reference_fasta = $options->{reference_fasta};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            my $interval_list = $options->{interval_list};
            my $bcftools = $options->{bcftools_exe};
            my $bcf_view_opts = $options->{bcftools_view_options};
	    my $min_recs = $self->options->{mimimum_calls};
	    
	    #my $max_cmdline_bams = $options->{max_cmdline_bams};
	    #if (scalar (@{$self->inputs->{bam_files}}) > $max_cmdline_bams) {
	    #	$self->warn("[todo] Generate a bam fofn");
	    #}
	    
	    my $req = $self->new_requirements(memory => 500, time => 1);
	    my $bam_list;
	    my ($bam_metadata,$basename);
	    
	    # if more than one bam, vcf basename and any meta data will be based upon the last one
            foreach my $bam (@{$self->inputs->{bam_files}}) {
		$bam_metadata = $bam->metadata;	
		$basename = $bam->basename;	
		my $bam_path = $bam->path;
		$bam_list .= "$bam_path ";
            }
	    
	    $mpileup_opts .= " -l $interval_list " if $interval_list;
	    $basename =~ s/\.bam/.mpileup.vcf.gz/;
	    
	    my $vcf_file = $self->output_file(output_key => 'vcf_files', basename => $basename, type => 'vcf');
	    my $vcf_path = $vcf_file->path;
	    if ($bam_metadata) {
		$vcf_file->add_metadata($bam_metadata);
	    }
	    
	    my $mpileup_cmd = qq[$samtools mpileup $mpileup_opts -f $reference_fasta $bam_list | $bcftools view $bcf_view_opts - | bgzip -c > $vcf_path];
	    
	    my $cmd = "use VRPipe::Steps::mpileup_vcf; VRPipe::Steps::mpileup_vcf->run_mpileup('$mpileup_cmd','$min_recs');";
	    $self->dispatch_vrpipecode($cmd, $req, {output_files => [$vcf_file]});
	};
    }
    method description {
        return "Run samtools mpileup and bcftools for one or more bams, generating one vcf without an intermediate bcf";
    }

    method outputs_definition {
        return { vcf_files => VRPipe::StepIODefinition->get(type => 'vcf', max_files => -1, description => 'a vcf file for each set of one or more input bams') };
    }

    method run_mpileup (ClassName|Object $self: Str $cmd_line, Int $min_recs) {
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my ($output_path) = $cmd_line =~ /> (\S+)$/;
        my $output_file = VRPipe::File->get(path => $output_path);
        $output_file->update_stats_from_disc;
        my $output_recs = $output_file->num_records;
	
        if ($output_recs < $min_recs) {
            $output_file->unlink;
	    $self->throw("Output VCF has $output_recs data lines, fewer than expected minimum $min_recs");
        }
        else {
            return 1;
        }
    }
}

1;
