=head1 NAME

VRPipe::Steps::bwa_aln_bam - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Steps::bwa_aln_bam with VRPipe::StepRole {
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file to map against'),
                 bwa_aln_options => VRPipe::StepOption->create(description => 'options to bwa aln, excluding the input fastq, reference, -b, -0/1/2 and -f options',
                                                            optional => 1),
                 bwa_exe => VRPipe::StepOption->create(description => 'path to your bwa executable',
                                                    optional => 1,
                                                    default_value => 'bwa') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam',
                                                            max_files => -1,
                                                            description => 'fastq files, which will be alnd independently',
                                                            metadata => {reads => 'total number of reads (sequences)'}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $bwa_exe = $options->{bwa_exe};
            my $bwa_opts = $options->{bwa_aln_options};
            if ($bwa_opts =~ /$ref|-f|-b|-[012]|aln/) {
                $self->throw("bwa_aln_options should not include the reference or -f,-b,-0,-1,-2 options, or the aln sub command");
            }
            my $cmd = $bwa_exe.' aln '.$bwa_opts;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'bwa', version => VRPipe::StepCmdSummary->determine_version($bwa_exe, '^Version: (.+)$'), summary => 'bwa aln '.$bwa_opts.' -f $sai_file -b $reference_fasta $bam_file'));
            
            my $req = $self->new_requirements(memory => 2900, time => 4);
            foreach my $bam (@{$self->inputs->{bam_files}}) {
                foreach my $mode (0, 1, 2) {
                    my $sai_file = $self->output_file(output_key => 'bwa_sai_files',
                                                    basename => $bam->basename.'.'.$mode.'.sai',
                                                    type => 'bin',
                                                    metadata => {source_bam => $bam->path->stringify,
                                                                 paired => $mode,
                                                                 reference_fasta => $ref->stringify,
                                                                 reads => $bam->metadata->{reads}});
                    
                    my $this_cmd = $cmd." -b -$mode -f ".$sai_file->path.' '.$ref.' '.$bam->path;
                    $self->dispatch([$this_cmd, $req, {output_files => [$sai_file]}]);
                }
            }
        };
    }
    method outputs_definition {
        return { bwa_sai_files => VRPipe::StepIODefinition->create(type => 'bin',
                                                                max_files => -1,
                                                                description => 'output files of independent bwa aln calls on each input bam',
                                                                metadata => {source_bam => 'the bam file that was input to bwa aln to generate this sai file',
                                                                             paired => 'what kind of reads were used from the source_bam: 0=single end reads; 1=forward reads; 2=reverse reads',
                                                                             reference_fasta => 'the reference fasta that reads were aligned to',
                                                                             reads => 'the number of reads in the source_bam'}) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Aligns the input bam(s) with bwa to the reference";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
