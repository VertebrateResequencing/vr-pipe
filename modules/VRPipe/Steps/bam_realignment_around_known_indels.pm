
=head1 NAME

VRPipe::Steps::bam_realignment_around_known_indels - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2013 Genome Research Limited.

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

# java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \
#  -jar /path/to/GenomeAnalysisTK.jar \
#  -I<lane-level.bam>   \
#  -R /path/to/reference.fasta \
#  -T IndelRealigner \
#  -targetIntervals<intervalListFromStep1Above.intervals>   \
#  -o<realigned.bam>   \
#  --compress 0 \
#  -known /path/to/Devine_Mills.indel_sites.vcf.gz \
#  -known /path/to/1000Genomes.indel_sites.vcf.gz \
#  -model KNOWNS_ONLY \
#  -LOD 0.4

class VRPipe::Steps::bam_realignment_around_known_indels extends VRPipe::Steps::gatk {
    around options_definition {
        return {
            %{ $self->$orig },
            known_indels_for_realignment => VRPipe::StepOption->create(description => 'the -known option(s) for GATK RealignerTargetCreator and IndelRealigner which define known indel sites. Could be --DBSNP and -B options for older versions of GATK.'),
            bam_realignment_options      => VRPipe::StepOption->create(description => 'command line options for GATK IndelRealigner; excludes -known options which are set by another StepOption', optional => 1, default_value => '-LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing'),
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more bam files'
            ),
            bai_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                max_files   => -1,
                description => 'index files for the input bam files'
            ),
            intervals_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'GATK intervals file for known indel sites'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $known_indels   = $options->{known_indels_for_realignment};
            my $intervals_file = $self->inputs->{intervals_file}->[0]->path;
            
            my $realign_opts = $options->{bam_realignment_options};
            if ($realign_opts =~ /$ref|known|IndelRealigner/) {
                $self->throw("gatk_realign_options should not include the reference, known indel files or IndelRealigner task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) ' . $realign_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_base       = $bam->basename;
                my $bam_meta       = $bam->metadata;
                my $realigned_base = $bam_base;
                $realigned_base =~ s/bam$/realign.bam/;
                my @outfiles;
                my $realigned_bam_file = $self->output_file(
                    output_key => 'realigned_bam_files',
                    basename   => $realigned_base,
                    type       => 'bam',
                    metadata   => $bam_meta
                );
                push @outfiles, $realigned_bam_file;
                
                unless ($realign_opts =~ m/--disable_bam_indexing/) {
                    my $index_base = $realigned_base;
                    $index_base =~ s/bam$/bai/;
                    push @outfiles,
                      $self->output_file(
                        output_key => 'realigned_bam_index_files',
                        basename   => $index_base,
                        type       => 'bin',
                        metadata   => $bam_meta
                      );
                }
                
                my $temp_dir = $options->{tmp_dir} || $realigned_bam_file->dir;
                
                my $this_cmd = $self->gatk_prefix($req->memory, $temp_dir) . qq[ -T IndelRealigner -R $ref -I ] . $bam->path . qq[ -o ] . $realigned_bam_file->path . qq[ -targetIntervals $intervals_file $known_indels $realign_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_realignment_around_known_indels', 'realign_and_check', [$this_cmd, $req, { output_files => \@outfiles }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            realigned_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a name-sorted uncompressed bam file with improved alignments near indels'
            ),
            realigned_bam_index_files => VRPipe::StepIODefinition->create(
                type        => 'bin',
                min_files   => 0,
                max_files   => -1,
                description => 'index file for the indel realigned bam file'
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Realigns reads around known indels to improve subsequent variant calling, producing a name-sorted uncompressed bam";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method realign_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /-I (\S+) -o (\S+)/;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->_filetype->check_records_vs_input($in_file, $cmd_line);
    }
}

1;
