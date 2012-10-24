
=head1 NAME

VRPipe::Steps::fastq_import - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::fastq_import with VRPipe::StepRole {
    use File::Fetch;
    use Net::FTP::Robust;
    
    method options_definition {
        return {};
    }
    
    method inputs_definition {
        return {
            fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                description => 'fastq files',
                max_files   => -1,
                metadata    => {
                    remote_path  => 'the complete remote location of the file',
                    expected_md5 => 'the md5 checksum the file is supposed to have',
                    optional     => ['remote_path', 'expected_md5']
                },
                check_existence => 0
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $fq_file (@{ $self->inputs->{fastq_files} }) {
                my $ifile = $fq_file->path;
                
                my $meta        = $fq_file->metadata;
                my $remote_path = $meta->{remote_path};
                if ($remote_path && $remote_path ne $ifile) {
                    # download to the path of our input file, which doesn't
                    # exist yet
                    my $ofile = $self->output_file(output_key => 'local_fastq_files', output_dir => $ifile->dir, basename => $ifile->basename, type => 'fq', metadata => $fq_file->metadata);
                    $self->dispatch_vrpipecode(qq[use VRPipe::Steps::fastq_import; VRPipe::Steps::fastq_import->download_fastq(source => q[$remote_path], dest => q[$ifile]);], $req, { output_files => [$ofile] });
                }
                else {
                    # symlink our existing input file to the pipeline output dir
                    # so that if this step is restarted, we won't delete our
                    # input file
                    my $ofile = $self->output_file(output_key => 'local_fastq_files', basename => $ifile->basename, type => 'fq', metadata => $fq_file->metadata);
                    $fq_file->symlink($ofile);
                }
            }
        };
    }
    
    method outputs_definition {
        return {
            local_fastq_files => VRPipe::StepIODefinition->create(
                type        => 'fq',
                description => 'a fastq file on a local disc',
                max_files   => -1
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "If fastq files in the datasource are on an external ftp site, downloads them to local disc";
    }
    
    method max_simultaneous {
        return 50;
    }
    
    method download_fastq (ClassName|Object $self: Str :$source!, Str|File :$dest!) {
        my $fq_file      = VRPipe::File->get(path => $dest);
        my $meta         = $fq_file->meta;
        my $expected_md5 = $meta->{expected_md5};
        my $local_dir    = $fq_file->dir;
        my $out_file     = $fq_file->path;
        
        $fq_file->disconnect;
        
        my $ff        = File::Fetch->new(uri => $source);
        my $scheme    = $ff->scheme;
        my $host      = $ff->host;
        my $path      = $ff->path;
        my $basename  = $ff->file;
        my $full_path = $path . $basename;
        
        if ($scheme eq 'ftp') {
            # use Net::FTP::Robust, since it's potentially better
            my $ftp = Net::FTP::Robust->new(Host => $host);
            
            $ftp->get($full_path, $local_dir);
            
            $fq_file->update_stats_from_disc;
            unless ($fq_file->s) {
                $self->throw("After 10 automated attempts, failed to download $source at all");
            }
            
            if ($expected_md5) {
                my $ok = $self->verify_md5($out_file, $expected_md5);
                
                unless ($ok) {
                    my $tries = 0;
                    while (!$ok) {
                        $fq_file->unlink;
                        
                        $fq_file->disconnect;
                        $ftp->get($full_path, $local_dir);
                        $fq_file->update_stats_from_disc;
                        $ok = $self->verify_md5($out_file, $expected_md5);
                        
                        $tries++;
                        last if $tries == 3;
                    }
                }
                
                unless ($ok) {
                    $fq_file->unlink;
                    $self->throw("Tried downloading $source 3 times, but the md5 never matched '$expected_md5'");
                }
            }
        }
        else {
            $self->throw("'$scheme' downloads are not yet supported");
        }
    }
}

1;
