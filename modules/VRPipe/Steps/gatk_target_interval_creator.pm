=head1 NAME

VRPipe::Steps::gatk_target_interval_creator - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

# java -Xmx1g -jar /path/to/GenomeAnalysisTK.jar \
#  -T RealignerTargetCreator \
#  -R /path/to/reference.fasta \
#  -o /path/to/output.intervals \
#  -known /path/to/Devine_Mills.indel_sites.vcf.gz \
#  -known /path/to/1000Genomes.indel_sites.vcf.gz

class VRPipe::Steps::gatk_target_interval_creator extends VRPipe::Steps::gatk {
    around options_definition {
        return { %{$self->$orig},
                 known_indels_for_realignment => VRPipe::StepOption->create(description => 'the -known option(s) for GATK RealignerTargetCreator and IndelRealigner which define known indel sites. Could be --DBSNP and -B options for older versions of GATK.'),
                 target_intervals_options => VRPipe::StepOption->create(description => 'command line options for GATK RealignerTargetCreator; excludes -known options which are set by another StepOption', optional => 1)};
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $intervals_opts = $options->{target_intervals_options};
            if ($intervals_opts =~ /$ref|known|RealignerTargetCreator/) {
                $self->throw("target_intervals_options should not include the reference, known file options or RealignerTargetCreator task command");
            }
            
            # Determine basename for intervals file
            my $known_indels = $options->{known_indels_for_realignment};
            my @knowns = $known_indels =~ m/-\S+ (\S+)/g;
            @knowns || $self->throw('No known indel sites supplied');
            my @known_files = map { Path::Class::File->new($_) } @knowns;
            my $basename = join '_', map { $_->basename } @known_files;
            $basename =~ s/\.vcf(\.gz)?$//g;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'GenomeAnalysisTK', 
                                   version => $self->gatk_version(),
                                   summary => 'java $jvm_args -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference_fasta -o $intervals_file -known $known_indels_file(s) '.$intervals_opts));
            
            my $intervals_file = $self->output_file(output_key => 'intervals_file',
                                              output_dir => $known_files[0]->dir->stringify, # arbitrarily place intervals file in the directory of the first listed file
                                              basename => qq[$basename.intervals],
                                              type => 'txt',
                                              metadata => { known_files => join(',', @knowns) });
            
            my $req = $self->new_requirements(memory => 4500, time => 1);
            my $jvm_args = $self->jvm_args($req->memory);
            
            my $this_cmd = $self->java_exe.qq[ $jvm_args -jar ].$self->jar.qq[ -T RealignerTargetCreator -R $ref -o ].$intervals_file->path.qq[ $known_indels $intervals_opts];
            $self->dispatch([$this_cmd, $req, {block_and_skip_if_ok => 1}]);
        };
    }
    method outputs_definition {
        return { intervals_file => VRPipe::StepIODefinition->create(type => 'txt', 
                                                                 description => 'GATK intervals file for known indel sites',
                                                                 metadata => { known_files => 'comma separated list of known indel file(s) used to create the intervals file' }) };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Creates target intervals file for known indel sites";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
