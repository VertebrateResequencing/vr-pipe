
=head1 NAME

VRPipe::Steps::mpileup_bcf_hapmap - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

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
use VRPipe::Parser;

class VRPipe::Steps::mpileup_bcf_hapmap extends VRPipe::Steps::mpileup_bcf {
    around options_definition {
        return { %{ $self->$orig },
                 samtools_mpileup_options => VRPipe::StepOption->create(description   => 'options for samtools mpileup, excluding -l and -f (-g is required)',
                                                                        optional      => 1,
                                                                        default_value => '-ugDI -d 1000 -C50') };
    }
    
    method inputs_definition {
        return { bam_files   => VRPipe::StepIODefinition->create(type => 'bam', max_files   => -1, description => 'bam files for bcf production'),
                 hapmap_file => VRPipe::StepIODefinition->create(type => 'txt', description => 'hapmap sites (-l positions) file') };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $samtools     = $options->{samtools_exe};
            my $mpileup_opts = $options->{samtools_mpileup_options};
            $self->throw("-g is required since we must make bcf files") unless $mpileup_opts =~ /g/;
            $self->throw("-l and -f options must not be used") if $mpileup_opts =~ /-l|-f/;
            
            my $hapmap_path = $self->inputs->{hapmap_file}->[0]->path;
            $mpileup_opts .= " -l $hapmap_path";
            my $req = $self->new_requirements(memory => 3900, time => 1);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                # work out the expected sample
                my $bam_path = $bam->path;
                my $meta     = $bam->metadata;
                my $sample   = $meta->{sample};
                unless ($sample) {
                    my $parser = VRPipe::Parser->create('bam', { file => $bam_path });
                    my %rg_info = $parser->readgroup_info();
                    my %samples;
                    while (my ($rg, $info) = each %rg_info) {
                        my $this_sample = $info->{SM} || next;
                        $samples{$this_sample} = 1;
                    }
                    my @samples = keys %samples;
                    if (@samples == 1) {
                        $sample = $samples[0];
                        $bam->add_metadata({ sample => $sample });
                    }
                }
                $self->throw("Unable to obtain sample name for bam file $bam_path") unless $sample;
                my $individual = $meta->{individual} || $sample;
                
                my $bcf_file = $self->output_file(output_key => 'bcf_files_with_metadata',
                                                  basename   => $bam_path->basename . '.bcf',
                                                  type       => 'bin',
                                                  metadata   => { sample => $sample, individual => $individual, source_bam => $bam_path->stringify });
                my $bcf_path = $bcf_file->path;
                my $cmd      = qq[$samtools mpileup $mpileup_opts -f $ref $bam_path > $bcf_path];
                $self->dispatch([$cmd, $req, { output_files => [$bcf_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { bcf_files_with_metadata => VRPipe::StepIODefinition->create(type        => 'bin',
                                                                             max_files   => -1,
                                                                             description => 'bcf file produced for bam file using samtools',
                                                                             metadata    => {
                                                                                           sample     => 'name of expected sample',
                                                                                           individual => 'name of expected individual',
                                                                                           source_bam => 'name of bam file used to produce bcf file' }) };
    }
    
    method description {
        return "Run samtools mpileup for each bam using the hapmap sites file provided to generate one bcf file with associated sample name metadata per bam";
    }
}

1;
