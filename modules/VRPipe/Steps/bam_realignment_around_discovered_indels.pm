
=head1 NAME

VRPipe::Steps::bam_realignment_around_discovered_indels - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bam_realignment_around_discovered_indels extends VRPipe::Steps::gatk {
    use VRPipe::Steps::bam_realignment_around_known_indels;
    
    around options_definition {
        return { %{ $self->$orig }, gatk_indelrealigner_options => VRPipe::StepOption->create(description => 'command line options for GATK IndelRealigner', optional => 1, default_value => '--disable_bam_indexing'), };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more bam files'
            ),
            intervals_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                max_files   => -1,
                description => 'GATK intervals file for each input bam file',
                metadata    => { source_bam => 'bam file this intervals file was created for' }
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $realign_opts = $options->{gatk_indelrealigner_options};
            if ($realign_opts =~ /$ref|IndelRealigner/) {
                $self->throw("gatk_indelrealigner_options should not include the reference or IndelRealigner task command");
            }
            
            $self->set_cmd_summary(
                VRPipe::StepCmdSummary->create(
                    exe     => 'GenomeAnalysisTK',
                    version => $self->gatk_version(),
                    summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file ' . $realign_opts
                )
            );
            
            my $req = $self->new_requirements(memory => 4500, time => 2);
            
            my @intervals_files = @{ $self->inputs->{intervals_file} };
            my @bam_files       = @{ $self->inputs->{bam_files} };
            $self->throw("number of bam files did not match number of intervals files") unless @bam_files == @intervals_files;
            my %bams = map { $_->resolve->path->stringify => $_ } @bam_files;
            
            foreach my $intervals_file (@{ $self->inputs->{intervals_file} }) {
                my $source_bam_path = $intervals_file->metadata->{source_bam};
                my $bam = VRPipe::File->get(path => $source_bam_path)->resolve;
                $self->throw("source bam '$source_bam_path' of intervals file was not an input bam") unless exists $bams{ $bam->path->stringify };
                my $bam_base       = $bam->basename;
                my $bam_meta       = $bam->metadata;
                my $realigned_base = $bam_base;
                $realigned_base =~ s/bam$/realign.bam/;
                my $realigned_bam_file = $self->output_file(
                    output_key => 'realigned_bam_files',
                    basename   => $realigned_base,
                    type       => 'bam',
                    metadata   => $bam_meta
                );
                
                my $temp_dir = $options->{tmp_dir} || $realigned_bam_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $intervals_file_path = $intervals_file->path;
                my $this_cmd            = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar . qq[ -T IndelRealigner -R $ref -I ] . $bam->path . qq[ -o ] . $realigned_bam_file->path . qq[ -targetIntervals $intervals_file_path $realign_opts];
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bam_realignment_around_known_indels', 'realign_and_check', [$this_cmd, $req, { output_files => [$realigned_bam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return {
            realigned_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a name-sorted bam file with improved alignments near indels'
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
        return 0;            # meaning unlimited
    }
}

1;
