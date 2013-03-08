
=head1 NAME

VRPipe::Steps::fake_fastq_metadata - a step

=head1 DESCRIPTION

similar to fastq_metadata step, but doesn't really do anything useful, since it
does not run fastqcheck: for testing purposes on all machines

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::fake_fastq_metadata with VRPipe::StepRole {
    method options_definition {
        return {};
    }
    
    method inputs_definition {
        return { fastq_files => VRPipe::StepIODefinition->create(type => 'fq', description => 'fastq files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            foreach my $fq_file (@{ $self->inputs->{fastq_files} }) {
                # our output file is our input file
                my $ifile = $fq_file->path;
                $self->output_file(output_key => 'fastq_files_with_metadata', output_dir => $ifile->dir, basename => $ifile->basename, type => 'fq');
                
                # we expect the datasource to fill in most if not all metadata,
                # but some things may need calculating
                my $meta = $fq_file->metadata;
                unless ($meta->{bases} && $meta->{reads} && $meta->{avg_read_length}) {
                    # run fastqcheck to generate stats
                    # record the fact that we (probably) made the fastqcheck
                    # file, even though it isn't in our outputs_definition *** or should we be able to have optional outputs?
                    my $fqc_file = $self->output_file(basename => $ifile->basename . '.fastqcheck', type => 'txt', temporary => 1);
                    my $ofile = $fqc_file->path;
                    unless ($fqc_file->s) {
                        my $req = $self->new_requirements(memory => 500, time => 1);
                        $self->dispatch_wrapped_cmd('VRPipe::Steps::fake_fastq_metadata', 'fake_stats_from_fastqcheck', ["fastqcheck $ifile > $ofile", $req, { output_files => [$fqc_file] }]);
                    }
                }
                unless ($fq_file->md5) {
                    $self->dispatch_md5sum($fq_file, $meta->{expected_md5});
                }
            }
        };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub {
            my $self = shift;
            
            my $all_ok = 1;
            foreach my $ofile (@{ $self->outputs->{fastq_files_with_metadata} }) {
                my $meta = $ofile->metadata;
                unless ($meta->{bases} && $meta->{reads} && $meta->{avg_read_length}) {
                    $self->warn("some metadata was missing from the dataelement of " . $ofile->path . ", and our fastqcheck run failed to resolved that");
                    $all_ok = 0;
                }
                
                unless ($ofile->md5) {
                    $self->warn("missing md5 for " . $ofile->path);
                    $all_ok = 0;
                }
            }
            
            return $all_ok;
        };
    }
    
    method description {
        return "Takes a fastq file and associates metadata with the file in the VRPipe database, making the fastq file usable in other fastq-related Steps";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method fake_stats_from_fastqcheck (ClassName|Object $self: Str $cmd_line) {
        my ($fq_path, $fqc_path) = $cmd_line =~ /^fastqcheck (\S+) > (\S+)$/;
        $fqc_path || $self->throw("bad cmd line [$cmd_line]");
        my $fq_file  = VRPipe::File->get(path => $fq_path);
        my $fqc_file = VRPipe::File->get(path => $fqc_path);
        
        if ($fq_path =~ /\.gz/) {
            $cmd_line = qq{sleep 5; gunzip -c $fq_path | cat > $fqc_path};
        }
        
        $fqc_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $fqc_file->update_stats_from_disc(retries => 3);
        if ($fqc_file->s) {
            my $new_meta;
            $new_meta->{bases}           = 500;
            $new_meta->{reads}           = 50;
            $new_meta->{avg_read_length} = 10;
            
            $fq_file->add_metadata($new_meta);
        }
        else {
            $self->throw("$fqc_file failed to be made");
        }
    }
}

1;
