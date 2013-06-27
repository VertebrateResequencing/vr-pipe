
=head1 NAME

VRPipe::Steps::sga_preprocess_many - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::sga_preprocess_many extends VRPipe::Steps::sga_preprocess {
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            
            my $sga_exe  = $options->{sga_exe};
            my $sga_opts = $options->{sga_preprocess_options};
            my $compress = $options->{sga_preprocess_compress_fastq};
            if ($sga_opts =~ /preprocess/) {
                $self->throw("sga_preprocess_options should not include the preprocess subcommand");
            }
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'sga', version => VRPipe::StepCmdSummary->determine_version($sga_exe, '^Version: (.+)$'), summary => 'sga preprocess ' . $sga_opts . ' $fastq_file(s)'));
            
            my %fastqs;
            foreach my $fq (@{ $self->inputs->{fastq_files} }) {
                my $meta = $fq->metadata;
                next unless $meta->{paired};
                my $basename = $fq->basename;
                $basename =~ s/_(1|2|M)\.(fq|fastq)(\.gz)?$/\.processed.fq/;
                if ($compress) {
                    $basename .= '.gz';
                }
                if ($meta->{paired} == 1) {
                    unshift @{ $fastqs{$basename} }, $fq;
                }
                else {
                    push @{ $fastqs{$basename} }, $fq;
                }
            }
            
            my $req = $self->new_requirements(memory => 3900, time => 1);
            my (@processed_fqs, @files);
            my $base_cmd = qq[$sga_exe preprocess $sga_opts ];
            foreach my $fq (keys %fastqs) {
                my @fq_ids = map { $_->id } @{ $fastqs{$fq} };
                my $meta = $self->common_metadata($fastqs{$fq});
                push @processed_fqs, $self->output_file(output_key => 'preprocessed_fastq_files', basename => $fq, type => 'fq', metadata => $meta);
                push @files, $processed_fqs[-1]->id . qq[=>[qw(@fq_ids)]];
            }
            my $files_string = join ',', @files;
            my $this_cmd = "use VRPipe::Steps::sga_preprocess_many; VRPipe::Steps::sga_preprocess_many->preprocess_many('$base_cmd', compress => $compress, files => { $files_string });";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@processed_fqs });
        };
    }
    
    method description {
        return "Prepare fastq files for assembly with sga.";
    }
    
    method preprocess_many (ClassName|Object $self: Str $base_cmd, Bool :$compress?, HashRef :$files) {
        keys %{$files} > 0 || $self->throw("You must supply file ids to this method");
        my $write_to = $compress ? '| gzip -c >' : '>';
        while (my ($output, $inputs) = each %{$files}) {
            my @fqs = map { VRPipe::File->get(id => $_)->path } @$inputs;
            my $output_path = VRPipe::File->get(id => $output)->path;
            my $cmd_line = $base_cmd . join(' ', @fqs) . " $write_to " . $output_path;
            my $ok = $self->preprocess_and_check($cmd_line);
            unless ($ok) {
                $self->throw("cmd [$cmd_line] failed");
            }
        }
        return 1;
    }
}

1;
