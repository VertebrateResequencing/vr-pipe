=head1 NAME

VRPipe::Steps::irods - a step

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

class VRPipe::Steps::irods with VRPipe::StepRole {
    has 'irods_exes' => (is => 'ro',
                         isa => 'HashRef',
                         lazy => 1,
                         builder => '_build_irods_exes');
    
    method _build_irods_exes {
        return {};
    }
    
    method handle_exes (HashRef $options) {
        my $irods_exes = $self->irods_exes;
        foreach my $exe (keys %{$irods_exes}) {
            my $method = $exe.'_exe';
            next unless defined $options->{$method};
            $irods_exes->{$exe} = $options->{$method};
        }
    }
    
    method options_definition {
        my %opts;
        
        my $irods_exes = $self->irods_exes;
        foreach my $exe (keys %{$irods_exes}) {
            $opts{$exe.'_exe'} = VRPipe::StepOption->get(description => "path to your irods '$exe' executable", optional => 1, default_value => $irods_exes->{$exe} || $exe);
        }
        
        return \%opts;
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub { return 1; };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generic step for steps using irods";
    }
    method max_simultaneous {
        return 15;
    }
    
    method get_file_by_basename (ClassName|Object $self: Str :$basename!, Str|File :$dest!, Str :$zone = 'seq',
                                 Str|File :$iget!, Str|File :$iquest!, Str|File :$ichksum!) {
        my $dest_file = VRPipe::File->get(path => $dest);
        $dest_file->disconnect;
        
        my $cmd = qq[$iquest -z $zone "SELECT COLL_NAME, DATA_NAME WHERE DATA_NAME = '$basename'"];
        open(my $irods, "$cmd |");
        my ($path, $filename);
        while (<$irods>) {
            chomp;
            # output looks like:
            # Zone is seq
            # COLL_NAME = /seq/5150
            # DATA_NAME = 5150_1#1.bam
            # -------------------------
            # 
            # or an error is thrown, and only:
            # Zone is seq
            # is output
            
            if (/^COLL_NAME = (.+)$/){
                $path = $1;
            }
            if (/^DATA_NAME = (.+)$/){
                $filename = $1;
            }
        }
        close $irods;
        
        if ($path && $filename) {
            unless ($filename eq $basename){
                $self->throw("$filename should be the same as $basename");
            }
            my $irods_file = join('/', ($path, $filename));
            
            $self->get_file(source => $irods_file, dest_file => $dest_file, iget => $iget, ichksum => $ichksum);
        }
        else {
            $self->throw("A file with basename $basename could not be found in iRods zone $zone");
        }
    }
    
    method get_file_md5 (ClassName|Object $self: Str :$file!, Str|File :$ichksum!) {
        my $md5 = `$ichksum $file`;
        chomp $md5;
        $md5 =~s/.*\s//;
        return $md5;
    }
    
    method get_file (ClassName|Object $self: Str :$source!, VRPipe::File :$dest_file!,
                     Str|File :$iget!, Str|File :$ichksum!) {
        my $dest = $dest_file->path;
        $dest_file->disconnect;

        # before we go fetch a file, check the md5 matches what we're expecting
        my $irodschksum = $self->get_file_md5(file => $source, ichksum => $ichksum);
        my $expected_md5 = $dest_file->metadata->{expected_md5} || $irodschksum;
        unless ($irodschksum eq $expected_md5) {
            $dest_file->unlink;
            $self->throw("expected md5 checksum in metadata did not match md5 of $source in IRODS; aborted");
        }
        
        # -K: checksum
        # -Q: use UDP rather than TCP
        # -f: force overwrite
        my $failed = system($iget, "-K", "-f", $source, $dest);
        
        $dest_file->update_stats_from_disc;
        if ($failed) {
            $dest_file->unlink;
            $self->throw("iget failed for $source -> $dest");
        }
        
        # correct permissions, just in case
        chmod 0664, $dest;
        
        # double-check the md5 (iget -K doesn't always work?)
        my $ok = $dest_file->verify_md5($dest, $expected_md5);
        unless ($ok) {
            $dest_file->unlink;
            $self->throw("we got $source -> $dest, but the md5 checksum did not match; deleted");
        }
    }
}

1;
