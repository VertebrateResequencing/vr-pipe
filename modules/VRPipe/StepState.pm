
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
PipelineSetup. It also provides the submissions() method to get at all the
submissions for a particular DataElement/PipelineSetup/Step combo.

Different DataElements for the same PipelineSetup, running the same Step, may
also end up wanting to run the exact same command line(s) (eg. when they're
running block_and_skip_if_ok jobs, like indexing a reference file in the first
step of a pipeline)). To avoid the wasteful creation of duplicate submissions
for the same job, same_submissions_as() can store another StepState that the
submissions were originally created for, and then submissions() will return
that StepState's submissions. submission_search_id() can be used when
constructing a Submission search query limited to the submissions of a
particular StepState, and it returns the appropriate StepState id().

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2014 Genome Research Limited.

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
    has 'stepmember' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::StepMember'
    );
    
    has 'dataelement' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::DataElement'
    );
    
    has 'pipelinesetup' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::PipelineSetup'
    );
    
    has 'cmd_summary' => (
        is          => 'rw',
        isa         => Persistent,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1,
        belongs_to  => 'VRPipe::StepCmdSummary'
    );
    
    has 'same_submissions_as' => (
        is          => 'rw',
        isa         => Persistent,
        coerce      => 1,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1,
        belongs_to  => 'VRPipe::StepState'
    );
    
    has 'complete' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    has 'dispatched' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 0
    );
    
    __PACKAGE__->make_persistent(has_many => [[submissions => 'VRPipe::Submission'], ['_output_files' => 'VRPipe::StepOutputFile']]);
    
    around same_submissions_as {
        shift;
        shift;
        if (@_) {
            my $id = shift;
            $self->throw("You can't set same_submissions_as to the same state as itself") if $id && $id == $self->id;
            return $self->$orig($id);
        }
        return $self->$orig;
    }
    
    around submissions {
        if (my $other_state = $self->same_submissions_as) {
            return $other_state->submissions;
        }
        return $self->$orig;
    }
    
    method submission_search_id {
        if (my $other_state = $self->same_submissions_as) {
            return $other_state->id;
        }
        return $self->id;
    }
    
    method output_files (Maybe[PersistentFileHashRef] $new_hash?, Bool :$only_unique_to_us?) {
        $only_unique_to_us = 0 if $new_hash;
        my @current_sofiles = VRPipe::StepOutputFile->search({ stepstate => $self->id, output_key => { '!=' => 'temp' } }, { prefetch => 'file' });
        my %hash;
        foreach my $sof (@current_sofiles) {
            my $file = $sof->file;
            
            if ($only_unique_to_us) {
                my $others = VRPipe::StepOutputFile->search({ file => $file->id, output_key => { '!=' => 'temp' }, stepstate => { '!=' => $self->id }, 'stepstate.pipelinesetup' => { '!=' => $self->pipelinesetup->id } }, { join => { stepstate => 'pipelinesetup' } });
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
                    $self->pipelinesetup->log_event("Deleting StepOutputFile for " . VRPipe::File->get(id => $file_id)->path . " because output_files() call with a new hash did not have this file in it", dataelement => $self->dataelement, stepstate => $self->id);
                    VRPipe::StepOutputFile->get(stepstate => $self, file => $file_id, output_key => $key)->delete;
                }
            }
            
            # remember new ones. *** there's some odd issue where the bulk
            # create fails because the constraint on file ids fails - because by
            # the time we bulk create the file rows have been changed?? To try
            # and avoid that we lock the files in a transaction
            delete $new_hash->{temp};
            my $transaction = sub {
                my @sof_args;
                my $own_id = $self->id;
                while (my ($key, $files) = each %$new_hash) {
                    foreach my $file (@$files) {
                        my ($current_file) = VRPipe::File->search({ id => $file->id }, { for => 'update' });
                        if (!$current_file || $current_file->path ne $file->path) {
                            $current_file = VRPipe::File->create(
                                path     => $file->path,
                                metadata => $file->metadata,
                                type     => $file->type,
                                e        => $file->e,
                                s        => $file->s
                            );
                            ($current_file) = VRPipe::File->search({ id => $current_file->id }, { for => 'update' });
                        }
                        
                        push(@sof_args, { stepstate => $own_id, file => $current_file->id, output_key => $key });
                    }
                }
                VRPipe::StepOutputFile->bulk_create_or_update(@sof_args);
            };
            $self->do_transaction($transaction, "StepOutputFile bulk_create_or_update failed");
            
            return $new_hash;
        }
        else {
            return \%hash;
        }
    }
    
    method temp_files (PersistentArrayRef $new_array?) {
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
    
    method start_over (Bool :$only_failed = 0) {
        my $log_self = 'StepState->start_over';
        $log_self .= '(only_failed => 1)' if $only_failed;
        $self->pipelinesetup->log_event("$log_self was called", stepstate => $self->id, dataelement => $self->dataelement->id, record_stack => 1);
        
        $self->block_until_locked;
        $self->maintain_lock;
        
        # before doing anything, record which output files we'll need to delete,
        # but don't actually delete anything until we've done all the db
        # updates
        my @files_to_unlink;
        unless ($only_failed) {
            foreach my $file ($self->output_files_list(only_unique_to_us => 1)) {
                push(@files_to_unlink, $file);
            }
        }
        
        # we need to be in a transaction or we risk leaving our selves in a
        # partial invalid state with eg. no submissions, no output files, but
        # claiming we are complete
        my $transaction = sub {
            # first remove output file rows from the db; we do this before
            # deleting subs to avoid a race condition where delete subs
            # while another process is parsing this, sees no subs and does a
            # parse and creates new sofs immediately before we then delete
            # them
            unless ($only_failed) {
                foreach my $sof ($self->_output_files) {
                    $self->pipelinesetup->log_event("$log_self call deleting StepOutputFile row for " . $sof->file->path, stepstate => $self->id, dataelement => $self->dataelement->id);
                    $sof->delete;
                }
            }
            
            # reset all associated submissions in order to reset their jobs
            foreach my $sub ($self->submissions) {
                if ($only_failed) {
                    next unless $sub->_failed();
                }
                
                $self->pipelinesetup->log_event("Calling Submission->start_over during a $log_self", stepstate => $self->id, dataelement => $self->dataelement->id, submission => $sub->id, job => $sub->job->id);
                $sub->start_over;
                
                # delete any stepstats there might be for us
                foreach my $ss (VRPipe::StepStats->search({ submission => $sub->id })) {
                    $ss->delete;
                }
                
                $self->pipelinesetup->log_event("$log_self call deleting Submission", stepstate => $self->id, dataelement => $self->dataelement->id, submission => $sub->id, job => $sub->job->id);
                $sub->delete;
            }
            
            # clear the dataelementstate to 0 steps completed; not important to try
            # and figure out the correct number of steps to set it to
            VRPipe::DataElementState->get(pipelinesetup => $self->pipelinesetup, dataelement => $self->dataelement, completed_steps => 0);
            
            # now reset self
            $self->complete(0);
            $self->dispatched(0);
            $self->update;
        };
        $self->do_transaction($transaction, "$log_self for " . $self->id . " failed");
        
        # now unlink the output files
        foreach my $file (@files_to_unlink) {
            $file->unlink;
        }
        
        $self->unlock;
        
        $self->pipelinesetup->log_event("$log_self call returning", stepstate => $self->id, dataelement => $self->dataelement->id);
    }
}

1;
