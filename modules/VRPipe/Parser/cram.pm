
=head1 NAME

VRPipe::Parser::cram - parse cram file headers

=head1 SYNOPSIS
    
    use VRPipe::Parser;
    
    # create object, supplying cram file
    my $pars = VRPipe::Parser->create('cram', {file => $cram_file});
    
    # get header information
    my $program = $pars->program();
    my %readgroup_info = $pars->readgroup_info();
    # etc.
    
=head1 DESCRIPTION

...

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Parser::cram with VRPipe::ParserRole {
    my $samtools_exe = file($ENV{SAMTOOLS}, 'samtools');
    
    has '_header' => (
        is      => 'rw',
        isa     => 'HashRef',
        default => sub { {} }
    );
    
    sub BUILD {
        my $self = shift;
        
        # ensure that _get_fh is called, which might be bypassed if user calls
        # no header methods
        $self->fh;
    }
    
    method _get_fh {
        my $file = $self->file;
        if (ref($file) && (File->check($file) || $file->isa('VRPipe::File'))) {
            my $vrpf = $file->isa('VRPipe::File') ? $file : VRPipe::File->create(path => $file->absolute, FileType->check($self->type) ? (type => $self->type) : ());
            $self->_set_vrpipe_file($vrpf); # just to keep hold of the filename
            my $filename = $vrpf->path->stringify;
            $self->{_filename} = $filename;
            
            # set up the open command which is just for the header
            my $open = "$samtools_exe view -H $filename |";
            open(my $fh, $open) || $self->throw("Couldn't open '$open': $!");
            
            $vrpf->disconnect;
            return $fh;
        }
        else {
            $self->throw("cram parser only supports VRPipe::File or path inputs");
        }
    }

=head2 sam_version
 
 Title   : sam_version
 Usage   : my $sam_version = $obj->sam_version();
 Function: Return the file format version of this sam file, as given in the
           header.
 Returns : number (undef if no header)
 Args    : n/a

=cut
    
    method sam_version {
        return $self->_get_single_header_tag('HD', 'VN');
    }

=head2 group_order
 
 Title   : group_order
 Usage   : my $group_order = $obj->group_order();
 Function: Return the group order of this sam file, as given in the header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut
    
    method group_order {
        return $self->_get_single_header_tag('HD', 'GO');
    }

=head2 sort_order
 
 Title   : sort_order
 Usage   : my $sort_order = $obj->sort_order();
 Function: Return the sort order of this sam file, as given in the header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut
    
    method sort_order {
        return $self->_get_single_header_tag('HD', 'SO');
    }

=head2 program_info
 
 Title   : program_info
 Usage   : my %all_program_info = $obj->program_info();
 Function: Get information about the programs used to create/process this bam,
           as reported in the header.
 Returns : undef if no PG lines in header, else:
           with no args: hash (keys are program ids, values are hash refs with
                               keys as tags (like VN and CL))
           with just a program id: hash (keys as tags, like VN and CL)
           with a program and a tag: the value of that tag for that program
 Args    : none for all info,
           program id for all the info for just that program,
           program id and tag (like 'VN' or 'CL') for specific info

=cut
    
    sub program_info {
        my ($self, @args) = @_;
        return $self->_handle_multi_line_header_types('PG', @args);
    }

=head2 program
 
 Title   : program
 Usage   : my $program = $obj->program();
 Function: Return the program used to do the mapping, as given in the header.
           If there is more than 1 PG header line, tries to guess which one is
           for the mapping program.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut
    
    method program {
        return $self->_guess_mapping_program();
    }
    
    method _guess_mapping_program {
        my %info     = $self->program_info();
        my @programs = keys %info;
        
        if (@programs == 1) {
            return $programs[0];
        }
        else {
            my (@known_prg, @unknown_prg);
            for my $program (@programs) {
                if ($program =~ /bwa|maq|ssaha|bfast|stampy|smalt/i) {
                    if (exists $info{$program}{PN}) {
                        push @known_prg, $info{$program}{PN};
                    }
                    else {
                        push @known_prg, $program;
                    }
                }
                elsif ($program !~ /GATK/) {
                    if (exists $info{$program}{PN}) {
                        push @unknown_prg, $info{$program}{PN};
                    }
                    else {
                        push @unknown_prg, $program;
                    }
                }
            }
            
            if (@known_prg) {
                return $known_prg[0];
            }
            elsif (@unknown_prg) {
                return $unknown_prg[0];
            }
            
            # guess randomly
            if (@programs) {
                return $programs[0];
            }
            else {
                # OMG, there's no PG lines in this bam file!
                return 'unknown_algorithm';
            }
        }
    }

=head2 program_version
 
 Title   : program_version
 Usage   : my $program_version = $obj->program_version();
 Function: Return the program version used to do the mapping, as given in the
           header.
           If there is more than 1 PG header line, tries to guess which one is
           for the mapping program.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut
    
    method program_version {
        my $program_id = $self->_guess_mapping_program();
        return $self->program_info($program_id, 'VN');
    }

=head2 command_line
 
 Title   : command_line
 Usage   : my $command_line = $obj->command_line();
 Function: Return the command line used to do the mapping, as given in the
           header.
           If there is more than 1 PG header line, tries to guess which one is
           for the mapping program.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut
    
    method command_line {
        my $program_id = $self->_guess_mapping_program();
        return $self->program_info($program_id, 'CL');
    }

=head2 sequence_info
 
 Title   : sequence_info
 Usage   : my %all_sequences_info = $obj->sequence_info();
           my %sequence_info = $obj->sequence_info('chr1');
           my $seq_length = $obj->sequence_info('chr1', 'LN');
 Function: Get information about the reference sequences, as reported in the
           header.
 Returns : undef if no SQ lines in header, else:
           with no args: hash (keys are sequence ids, values are hash refs with
                               keys as tags (like LN and M5))
           with just a sequence id: hash (keys as tags, like LN and M5)
           with a sequence and a tag: the value of that tag for that sequence
 Args    : none for all info,
           sequence id for all the info for just that sequence,
           sequence id and tag (like 'LN' or 'M5') for specific info

=cut
    
    sub sequence_info {
        my ($self, @args) = @_;
        return $self->_handle_multi_line_header_types('SQ', @args);
    }

=head2 readgroup_info
 
 Title   : readgroup_info
 Usage   : my %all_rg_info = $obj->readgroup_info();
           my %rg_info = $obj->readgroup_info('SRR00001');
           my $library = $obj->readgroup_info('SRR00001', 'LB');
 Function: Get information about the read groups, as reported in the header.
 Returns : undef if no RG lines in header, else:
           with no args: hash (keys are readgroups, values are hash refs with
                               keys as tags (like LB and SM))
           with just a readgroup id: hash (keys as tags, like LB and SM)
           with a readgroup and a tag: the value of that tag for that readgroup
 Args    : none for all info,
           readgroup id for all the info for just that readgroup,
           readgroup id and tag (like 'LB' or 'SM') for specific info

=cut
    
    sub readgroup_info {
        my ($self, @args) = @_;
        return $self->_handle_multi_line_header_types('RG', @args);
    }

=head2 samples
 
 Title   : samples
 Usage   : my @samples = $obj->samples();
 Function: Get all the unique SM fields from amongst all RG lines in
           the header.
 Returns : list of strings (sample names)
 Args    : none

=cut
    
    method samples {
        return $self->_get_unique_rg_fields('SM');
    }
    
    method _get_unique_rg_fields (Str $field) {
        my %vals;
        my %rg_info = $self->readgroup_info();
        while (my ($rg, $data) = each %rg_info) {
            $vals{ $data->{$field} || next } = 1;
        }
        my @uniques = sort keys %vals;
        return @uniques;
    }
    
    method _handle_multi_line_header_types (Str $type, Str $id?, Str $tag?) {
        my $lines = $self->_get_header_type($type) || return;
        
        # organise the data into by-id hash
        my %all_info;
        foreach my $line (@{$lines}) {
            my %this_data = $self->_tags_to_hash(@{$line});
            my $this_id = $this_data{SN} || $this_data{ID};
            delete $this_data{SN};
            delete $this_data{ID};
            
            $all_info{$this_id} = \%this_data;
        }
        
        if (defined $id) {
            my $id_info = $all_info{$id} || return;
            if ($tag) {
                return $id_info->{$tag};
            }
            else {
                return %{$id_info};
            }
        }
        else {
            return %all_info;
        }
    }
    
    method _get_single_header_tag (Str $type, Str $tag) {
        my $type_data = $self->_get_header_type($type) || return;
        
        my %data = $self->_tags_to_hash(@{$type_data});
        
        return $data{$tag};
    }
    
    method _tags_to_hash (@tags) {
        my %hash;
        foreach my $tag (@tags) {
            my ($this_tag, $value) = $tag =~ /^(\w\w):(.+)/;
            $hash{$this_tag} = $value;
        }
        return %hash;
    }
    
    method _get_header_type (Str $type) {
        my $fh = $self->fh() || return;
        
        $self->_get_header();
        my $hash = $self->_header;
        
        if (defined $hash->{$type}) {
            return $hash->{$type};
        }
        
        return;
    }
    
    method _get_header {
        my $fh = $self->fh() || return;
        return 1 if $self->_header_parsed();
        
        my $hash = $self->_header;
        
        while (<$fh>) {
            if (/^@/) {
                my @tags = split("\t", $_);
                my $type = shift @tags;
                $type = substr($type, 1);
                
                if ($type eq 'HD') {
                    # we only expect and handle one of these lines per file
                    $hash->{$type} = \@tags;
                }
                else {
                    push(@{ $hash->{$type} }, \@tags);
                }
            }
            else {
                $self->throw("Got a header line from a bam file that did not start with '\@'!: $_");
            }
        }
        
        $self->_header($hash);
        $self->_set_header_parsed();
        return 1;
    }
    
    sub close {
        my $self = shift;
        
        # make sure we've finished reading the whole thing before attempting to
        # close
        my $fh = $self->fh || return;
        while (<$fh>) {
            next;
        }
        close($fh);
    }
    
    method next_record {
        return 0;
    }
}

1;
