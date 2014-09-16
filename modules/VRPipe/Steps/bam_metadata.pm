
=head1 NAME

VRPipe::Steps::bam_metadata - a step

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

class VRPipe::Steps::bam_metadata extends VRPipe::Steps::bamcheck {
    has 'meta_to_check' => (
        is      => 'rw',
        isa     => 'ArrayRef',
        builder => '_build_meta_to_check'
    );
    
    my $core_samtools_exe = file($ENV{SAMTOOLS}, 'samtools');
    
    around options_definition {
        my $options = $self->$orig;
        # we ignore bamcheck options, since bam_metadata is used to get the
        # full basic stats on bams, regardless of if they are exome etc.
        delete $options->{bamcheck_options};
        delete $options->{reference_fasta};
        delete $options->{exome_targets_file};
        
        # and we add an option to avoid adding original_pg_chain, which is
        # inappropriate and possibly too large to store on eg. merged bams
        $options->{store_original_pg_chain} = VRPipe::StepOption->create(
            description   => 'If your input bam was not created by VRPipe and will subsequently go through the bam_reheader step, keep this on; otherwise be sure to turn it off.',
            optional      => 1,
            default_value => 1
        );
        
        $options->{bam_metadata_bamcheck_options} = VRPipe::StepOption->create(
            description => 'options to bam_metadata which will be used as arguments to bamcheck',
            optional    => 1
        );
        
        return $options;
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options       = $self->options;
            my $bamcheck_exe  = $options->{bamcheck_exe};
            my $store_pg      = $options->{store_original_pg_chain};
            my $bamcheck_opts = $options->{bam_metadata_bamcheck_options} ? $options->{bam_metadata_bamcheck_options} : '';
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $ifile = $bam_file->path;
                
                # run bamcheck if we don't have enough metadata
                my $meta       = $bam_file->metadata;
                my $meta_count = 0;
                foreach my $type (@{ $self->meta_to_check }) {
                    $meta_count++ if $meta->{$type};
                }
                unless ($meta_count == @{ $self->meta_to_check }) {
                    my $check_file = $self->output_file(basename => $ifile->basename . '.bamcheck', type => 'txt', temporary => 1);
                    my $ofile = $check_file->path;
                    $self->dispatch_wrapped_cmd('VRPipe::Steps::bamcheck', 'stats_from_bamcheck', ["$bamcheck_exe $bamcheck_opts $ifile > $ofile", $req, { output_files => [$check_file] }]);
                }
                
                # we'll also check the header for existing PG lines and store
                # those as metadata
                if ($store_pg && !defined $meta->{original_pg_chain}) {
                    $self->dispatch_vrpipecode("use VRPipe::Steps::bam_metadata; VRPipe::Steps::bam_metadata->store_pg_chain(bam => q[$ifile]);", $req);
                }
            }
        };
    }
    
    method _build_meta_to_check {
        my @meta_to_check = (qw(bases reads avg_read_length forward_reads reverse_reads rmdup_reads));
        return \@meta_to_check;
    }
    
    method outputs_definition {
        return {};
    }
    
    method description {
        return "Takes a bam file and associates metadata with the file in the VRPipe database, making the bam file usable in other bam-related Steps";
    }
    
    method store_pg_chain (ClassName|Object $self: Str|File :$bam!) {
        my $open = "$core_samtools_exe view -H $bam |";
        open(my $fh, $open) || $self->throw("Couldn't open '$open': $!");
        my $pg_chain;
        while (<$fh>) {
            next unless /^\@PG/;
            $pg_chain .= $_;
        }
        close($fh);
        
        if ($pg_chain && (length($pg_chain) < 20000)) {
            my $bam_file = VRPipe::File->get(path => $bam);
            $bam_file->add_metadata({ original_pg_chain => $pg_chain });
        }
        return 1;
    }
}

1;
