
=head1 NAME

VRPipe::Steps::bam_to_fastq_many - a step

=head1 DESCRIPTION

Runs the bam2fastq executable to extract sequences from a BAM file into fastq

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Steps::bam_to_fastq_many extends VRPipe::Steps::bam_to_fastq {
    method body_sub {
        return sub {
            my $self           = shift;
            my $options        = $self->options;
            my $bam2fastq_exe  = $options->{bam2fastq_exe};
            my $bam2fastq_opts = $options->{bam2fastq_opts};
            my $fastqcheck_exe = $options->{fastqcheck_exe};
            $bam2fastq_opts .= ' -f '; # force overwrite
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            
            my (@outfiles, @io_map);
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $meta = $bam->metadata;
                next unless $meta->{reads};
                my $paired = $meta->{paired};
                
                my $source_bam = $bam->path->stringify;
                my $fastq_meta = { source_bam => $source_bam };
                foreach my $key (qw(lane insert_size mean_insert_size library sample center_name platform study split_sequence population continent chrom from to individual project species)) {
                    if (defined $meta->{$key}) {
                        $fastq_meta->{$key} = $meta->{$key};
                    }
                }
                my $basename = $bam->basename;
                $basename =~ s/\.bam//;
                # '#' and '%' are special chars used by fastq output file definition
                # so must be replaced when generating fastq filenames from the bam filenames
                $basename =~ s/[#%]/_/g;
                my $out_spec;
                
                if ($paired) {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "${basename}_1.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads  => $meta->{forward_reads},
                            paired => 1
                        }
                    );
                    my $reverse = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "${basename}_2.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads  => $meta->{reverse_reads},
                            paired => 2
                        }
                    );
                    push @outfiles, ($fastq, $reverse);
                    
                    $fastq->add_metadata({ mate => $reverse->path->stringify });
                    $reverse->add_metadata({ mate => $fastq->path->stringify });
                    
                    push @io_map, $bam->id . '=>[' . $fastq->id . ',' . $reverse->id . ']';
                }
                else {
                    my $fastq = $self->output_file(
                        output_key => 'fastq_files',
                        basename   => "${basename}_M.fastq",
                        type       => 'fq',
                        metadata   => {
                            %$fastq_meta,
                            reads           => $meta->{reads},
                            bases           => $meta->{bases},
                            avg_read_length => sprintf("%0.2f", $meta->{reads} / $meta->{bases}),
                            paired          => 0
                        }
                    );
                    push @outfiles, $fastq;
                    
                    push @io_map, $bam->id . '=>[' . $fastq->id . ']';
                }
                
                my $out_log = $self->output_file(temporary => 1, basename => "$basename.log", type => 'txt');
                push(@outfiles, $out_log);
            }
            
            my $files_string = join ',', @io_map;
            my $this_cmd = "use VRPipe::Steps::bam_to_fastq_many; VRPipe::Steps::bam_to_fastq_many->bam_to_fastq_many(files => { $files_string }, bam2fastq_exe => q[$bam2fastq_exe], bam2fastq_opts => q[$bam2fastq_opts], fastqcheck_exe => q[$fastqcheck_exe]);";
            $self->dispatch_vrpipecode($this_cmd, $req, { output_files => \@outfiles });
        };
    }
    
    method bam_to_fastq_many (ClassName|Object $self: HashRef :$files!, Str|File :$bam2fastq_exe, Str :$bam2fastq_opts, Str|File :$fastqcheck_exe) {
        keys %{$files} > 0 || $self->throw("You must supply file ids to this method");
        while (my ($input, $outputs) = each %{$files}) {
            my @fqs = map { VRPipe::File->get(id => $_)->path } @$outputs;
            my $bam = VRPipe::File->get(id => $input)->path;
            my %output_hash;
            if (scalar @fqs == 1) {
                $output_hash{single} = $fqs[0];
            }
            elsif (scalar @fqs == 2) {
                $output_hash{forward} = $fqs[0];
                $output_hash{reverse} = $fqs[1];
            }
            else {
                $self->throw("Only 1 or 2 fastq files should be supplied");
            }
            my $ok = $self->bam_to_fastq(bam => $bam, %output_hash, bam2fastq_exe => $bam2fastq_exe, bam2fastq_opts => $bam2fastq_opts, fastqcheck_exe => $fastqcheck_exe);
            unless ($ok) {
                $self->throw("bam_to_fastq failed for $bam");
            }
        }
        return 1;
    }
}

1;
