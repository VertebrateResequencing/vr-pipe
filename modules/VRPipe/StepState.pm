
=head1 NAME

VRPipe::StepState - track Step completion for a particualar setup and element

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

Multiple different L<VRPipe::PipelineSetup>s can use the same
L<VRPipe::DataSource>, and even the same L<VRPipe::Pipeline>, and so the same
L<VRPipe::Step>s could be processing the same input files specified by the same
L<VRPipe::DataElement>s. StepState lets B<VRPipe> keep track of when a Step
successfully completes its work on a particular DataElement for a particular
PipelineSetup.

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

class VRPipe::StepState extends VRPipe::Persistent {
    has 'stepmember' => (is         => 'rw',
                         isa        => Persistent,
                         coerce     => 1,
                         traits     => ['VRPipe::Persistent::Attributes'],
                         is_key     => 1,
                         belongs_to => 'VRPipe::StepMember');
    
    has 'dataelement' => (is         => 'rw',
                          isa        => Persistent,
                          coerce     => 1,
                          traits     => ['VRPipe::Persistent::Attributes'],
                          is_key     => 1,
                          belongs_to => 'VRPipe::DataElement');
    
    has 'pipelinesetup' => (is         => 'rw',
                            isa        => Persistent,
                            coerce     => 1,
                            traits     => ['VRPipe::Persistent::Attributes'],
                            is_key     => 1,
                            belongs_to => 'VRPipe::PipelineSetup');
    
    has 'cmd_summary' => (is          => 'rw',
                          isa         => Persistent,
                          coerce      => 1,
                          traits      => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1,
                          belongs_to  => 'VRPipe::StepCmdSummary');
    
    has 'complete' => (is      => 'rw',
                       isa     => 'Bool',
                       traits  => ['VRPipe::Persistent::Attributes'],
                       default => 0);
    
    __PACKAGE__->make_persistent(has_many => [[submissions => 'VRPipe::Submission'], ['_output_files' => 'VRPipe::StepOutputFile']]);
    
    method output_files (Maybe[PersistentFileHashRef] $new_hash?, Bool :$only_unique_to_us?) {
        $only_unique_to_us = 0 if $new_hash;
        my @current_sofiles = VRPipe::StepOutputFile->search({ stepstate => $self->id, output_key => { '!=' => 'temp' } }, { prefetch => 'file' });
        my %hash;
        foreach my $sof (@current_sofiles) {
            my $file = $sof->file;
            
            if ($only_unique_to_us) {
                my $others = VRPipe::StepOutputFile->search({ file => $file->id, output_key => { '!=' => 'temp' }, stepstate => { '!=' => $self->id } });
                next if $others;
            }
            
            push(@{ $hash{ $sof->output_key } }, $file);
        }
        
        if ($new_hash) {
            # forget output files we no longer have
            while (my ($key, $files) = each %hash) {
                my @current_file_ids = map { $_->id } @$files;
                my @files_to_forget;
                if (exists $new_hash->{$key}) {
                    my %new_file_ids = map { $_->id => 1 } @{ $new_hash->{$key} };
                    foreach my $id (@current_file_ids) {
                        unless (exists $new_file_ids{$id}) {
                            push(@files_to_forget, $id);
                        }
                    }
                }
                else {
                    @files_to_forget = @current_file_ids;
                }
                
                foreach my $file_id (@files_to_forget) {
                    VRPipe::StepOutputFile->get(stepstate => $self, file => $file_id, output_key => $key)->delete;
                }
            }
            
            # remember new ones
            delete $new_hash->{temp};
            my @sof_args;
            my $own_id = $self->id;
            while (my ($key, $files) = each %$new_hash) {
                foreach my $file (@$files) {
                    push(@sof_args, { stepstate => $own_id, file => $file->id, output_key => $key });
                }
            }
            VRPipe::StepOutputFile->bulk_create_or_update(@sof_args);
            
            return $new_hash;
        }
        else {
            return \%hash;
        }
    }
    
    method temp_files (ArrayRefOfPersistent $new_array?) {
        my @file_ids = VRPipe::StepOutputFile->get_column_values('file', { stepstate => $self->id, output_key => 'temp' }); # *** can we get file objects out efficiently?
        
        if ($new_array) {
            # forget temp files we no longer have
            my %new_file_ids = map { $_->id => 1 } @$new_array;
            foreach my $file_id (@file_ids) {
                unless (exists $new_file_ids{$file_id}) {
                    VRPipe::StepOutputFile->get(stepstate => $self, file => $file_id, output_key => 'temp')->delete;
                }
            }
            
            # remember new ones
            my @sof_args;
            my $own_id = $self->id;
            foreach my $file (@$new_array) {
                push(@sof_args, { stepstate => $own_id, file => $file->id, output_key => 'temp' });
            }
            VRPipe::StepOutputFile->bulk_create_or_update(@sof_args);
            
            return $new_array;
        }
        else {
            return [map { VRPipe::File->get(id => $_) } @file_ids];
        }
    
    }
    
    method output_files_list (Bool :$only_unique_to_us = 0) {
        my $outputs = $self->output_files(undef, only_unique_to_us => $only_unique_to_us);
        my @files;
        if ($outputs) {
            foreach my $val (values %$outputs) {
                push(@files, @$val);
            }
        }
        return @files;
    }
    
    method update_output_file_stats {
        foreach my $file ($self->output_files_list) {
            $file->update_stats_from_disc(retries => 3);
        }
    }
    
    method unlink_output_files (Bool :$only_unique_to_us = 1) {
        foreach my $file ($self->output_files_list(only_unique_to_us => $only_unique_to_us)) {
            $file->unlink;
        }
    }
    
    method unlink_temp_files {
        foreach my $vrfile (@{ $self->temp_files }) {
            $vrfile->unlink;
        }
    }
    
    method start_over {
        $self->debug("start_over called for stepstate " . $self->id);
        
        # first reset all associated submissions in order to reset their jobs
        foreach my $sub ($self->submissions) {
            $sub->start_over;
            
            # delete any stepstats there might be for us
            foreach my $ss (VRPipe::StepStats->search({ submission => $sub->id })) {
                $ss->delete;
            }
            
            $sub->delete;
        }
        
        # clear the dataelementstate to 0 steps completed; not important to try
        # and figure out the correct number of steps to set it to
        VRPipe::DataElementState->get(pipelinesetup => $self->pipelinesetup, dataelement => $self->dataelement, completed_steps => 0);
        
        # now reset self
        $self->unlink_output_files(only_unique_to_us => 1);
        $self->complete(0);
        $self->update;
        
        # remove output file rows from the db
        foreach my $sof ($self->_output_files) {
            $sof->delete;
        }
    }
}

1;
