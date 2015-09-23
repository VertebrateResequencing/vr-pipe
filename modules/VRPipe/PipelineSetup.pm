
=head1 NAME

VRPipe::PipelineSetup - set up a pipeline

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The main thing that users want to do is "run a pipeline" on a given set of
data. Pipelines are modelled by L<VRPipe::Pipline> objects, and the set of data
is modelled by a L<VRPipe::DataSource>. A PipelineSetup relates the two and
also stores the user configuration of both, resulting in a named object
defining what the user wants to happen.

So users "set up" pipelines and get back a PipelineSetup. They, and the
B<VRPipe> system itself, look to these PipelineSetups to see what is supposed
to be run and how, on what data.

Multiple PipelineSetups can run at once, even working on the same Pipelines and
DataSources (presumably with different configuration options), and the system
ensures that there are no problems with similar work being done by different
PipelineSetups overwriting each other's files.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2015 Genome Research Limited.

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

class VRPipe::PipelineSetup extends VRPipe::Persistent {
    use POSIX qw(ceil);
    use DateTime;
    use DateTime::Format::Natural;
    use DateTime::TimeZone;
    use VRPipe::Config;
    use VRPipe::MessageTracker;
    use Sys::Hostname;
    
    our $local_timezone = DateTime::TimeZone->new(name => 'local');
    our $log_host = hostname();
    
    has 'name' => (
        is     => 'rw',
        isa    => Varchar [128],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'datasource' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::DataSource'
    );
    
    has 'pipeline' => (
        is         => 'rw',
        isa        => Persistent,
        coerce     => 1,
        traits     => ['VRPipe::Persistent::Attributes'],
        is_key     => 1,
        belongs_to => 'VRPipe::Pipeline'
    );
    
    has 'output_root' => (
        is     => 'rw',
        isa    => Dir,
        coerce => 1,
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'options' => (
        is                   => 'rw',
        isa                  => 'HashRef',
        traits               => ['VRPipe::Persistent::Attributes'],
        default              => sub { {} },
        allow_key_to_default => 1,
        is_key               => 1
    );
    
    has 'description' => (
        is          => 'rw',
        isa         => Text,
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'active' => (
        is      => 'rw',
        isa     => 'Bool',
        traits  => ['VRPipe::Persistent::Attributes'],
        default => 1
    );
    
    has 'user' => (
        is      => 'rw',
        isa     => Varchar [64],
        traits  => ['VRPipe::Persistent::Attributes'],
        default => scalar(getpwuid($<))
    );
    
    has 'unix_group' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'desired_farm' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has 'controlling_farm' => (
        is          => 'rw',
        isa         => Varchar [64],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    __PACKAGE__->make_persistent(has_many => [states => 'VRPipe::StepState']);
    
    # because lots of frontends get the dataelementstates for a particular
    # pipelinesetup
    sub _des_search_args {
        my $self              = shift;
        my $include_withdrawn = shift;
        my $only_withdrawn    = shift;
        my @terms             = ();
        if ($only_withdrawn) {
            @terms = ('dataelement.withdrawn' => 1);
        }
        elsif (!$include_withdrawn) {
            @terms = ('dataelement.withdrawn' => 0);
        }
        return ({ pipelinesetup => $self->id, @terms }, { prefetch => 'dataelement' });
    }
    
    method dataelementstates_pager (Bool :$include_withdrawn = 0, Bool :$only_withdrawn = 0) {
        return VRPipe::DataElementState->search_paged($self->_des_search_args($include_withdrawn, $only_withdrawn));
    }
    
    method dataelementstates (Bool :$include_withdrawn = 0, Bool :$only_withdrawn = 0) {
        return VRPipe::DataElementState->search($self->_des_search_args($include_withdrawn, $only_withdrawn));
    }
    
    around desired_farm (Maybe[Str] $farm?) {
        my $current_farm = $self->$orig();
        if ($farm) {
            if ($current_farm && $farm ne $current_farm) {
                $self->controlling_farm(undef);
            }
            return $self->$orig($farm);
        }
        return $current_farm;
    }
    
    method currently_complete {
        # we consider ourselves incomplete if we have no dataelementstates
        # assigned to us
        my $num_elements = VRPipe::DataElementState->search({ pipelinesetup => $self->id });
        return 0 unless $num_elements;
        
        # and we're incomplete if any of them have not completed all steps in
        # our pipeline
        my $num_steps_to_complete = $self->pipeline->step_members;
        my $elements_incomplete = VRPipe::DataElementState->search({ pipelinesetup => $self->id, completed_steps => { '<' => $num_steps_to_complete }, 'dataelement.withdrawn' => 0 }, { join => 'dataelement' });
        return $elements_incomplete ? 0 : 1;
    }
    
    method trigger (Bool :$first_step_only = 0, Bool :$prepare_elements = 1, VRPipe::DataElement :$dataelement?, Bool :$debug = 0) {
        my $setup_id     = $self->id;
        my $pipeline     = $self->pipeline;
        my @step_members = $pipeline->step_members;
        my $num_steps    = scalar(@step_members);
        my $datasource   = $self->datasource;
        
        warn "in trigger\n" if $debug;
        
        my $output_root = $self->output_root;
        $self->make_path($output_root);
        
        # work out which step inputs each step needs
        my %step_inputs;
        foreach my $adaptor (@{ $pipeline->adaptors || [] }) {
            my $hash = $adaptor->adaptor_hash || next;
            my %from_steps;
            while (my ($to_key, $from_keys) = each %{$hash}) {
                foreach my $from_step (values %{$from_keys}) {
                    next unless $from_step;
                    $from_steps{$from_step} = 1;
                }
            }
            if (keys %from_steps) {
                $step_inputs{ $adaptor->to_step } = [sort { $a <=> $b } keys %from_steps];
            }
        }
        warn "worked out step inputs\n" if $debug;
        
        # we either loop through all incomplete elementstates, or the
        # (single) elementstate for the supplied dataelement
        my $pager;
        my $mode;
        if ($dataelement) {
            $mode = 'single de';
            $pager = VRPipe::DataElementState->search_paged({ pipelinesetup => $setup_id, dataelement => $dataelement->id, completed_steps => $first_step_only ? 0 : { '<', $num_steps }, 'dataelement.withdrawn' => 0 }, { prefetch => 'dataelement' });
        }
        else {
            $mode = 'all des';
            # incomplete_element_states will block_until_locked to do datasource
            # updates
            eval { $pager = $datasource->incomplete_element_states($self, prepare => $prepare_elements, only_not_started => $first_step_only, debug => $debug); };
            if ($@) {
                return "DataSource error: $@";
            }
        }
        warn "got ", $pager->total_entries, " incomplete element states in mode '$mode'\n" if $debug;
        return unless $pager; #*** is it ever an error to have no pager?
        
        my $error_message;
        my $max_redos = $num_steps;
        my %cached_sameas;
        while (my $estates = $pager->next) {
            foreach my $estate (@$estates) {
                my $element         = $estate->dataelement;
                my $completed_steps = $estate->completed_steps;
                my $esid            = $estate->id;
                warn " working on estate $esid, which has completed $completed_steps / $num_steps\n" if $debug;
                next if $completed_steps == $num_steps;
                
                # we use a redis lock to skip the whole des to avoid triggering
                # the same thing multiple times simultaneously. Since we just
                # want to skip and do nothing (even in the case of 'single de'
                # mode) if another process is looking at this des, we don't have
                # to do anything fancy
                unless ($debug || $estate->lock) {
                    warn " estate $esid is being triggered by another process, will skip it\n";
                    next;
                }
                
                my $sm_error;
                SSTATE: foreach my $member (@step_members) {
                    my $step_number = $member->step_number;
                    next unless $step_number > $estate->completed_steps;
                    
                    my $state = VRPipe::StepState->create(
                        stepmember    => $member,
                        dataelement   => $element,
                        pipelinesetup => $self
                    );
                    
                    warn "  working on stepstate ", $state->id, "\n" if $debug;
                    
                    my ($do_next, $do_last);
                    my %previous_step_outputs = ();
                    my $pso_calculated        = 0;
                    
                    for (1 .. $max_redos) {
                        my $transaction = sub {
                            unless ($pso_calculated) {
                                # work out the previous step outputs relevant to
                                # this step
                                foreach my $input_step_number (@{ $step_inputs{$step_number} || [] }) {
                                    my $input_member = $step_members[$input_step_number - 1];
                                    my ($input_ss) = VRPipe::StepState->search({ stepmember => $input_member, dataelement => $element, pipelinesetup => $self });
                                    unless ($input_ss) {
                                        die "DataElement ", $element->id, " passed step $input_step_number, but the StepState corresponding to that doesn't exist! You'll probably have to 'vrpipe-elements --setup ", $self->id, " --start_from_scratch ...' this and other stuck dataelements to redo step $input_step_number\n";
                                    }
                                    
                                    my $input_step = $input_member->step(step_state => $input_ss);
                                    while (my ($key, $val) = each %{ $input_step->outputs() }) {
                                        $previous_step_outputs{$key}->{$input_step_number} = $val;
                                    }
                                }
                                $pso_calculated = 1;
                            }
                            warn "   worked out previous_step_outputs\n" if $debug;
                            
                            my $step = $member->step(previous_step_outputs => \%previous_step_outputs, step_state => $state, debug => $debug);
                            if ($state->complete || ($state->same_submissions_as && $state->same_submissions_as->complete)) {
                                warn "   stepstate is complete\n" if $debug;
                                $self->_complete_state($step, $state, $step_number, $pipeline, $estate);
                                $do_next = 1;
                                return;
                            }
                            
                            undef $sm_error;
                            my $error_ident = 'step ' . $step->name . " for setup $setup_id (dataelement " . $element->id . ', stepstate ' . $state->id . ')';
                            
                            # have we previously done the dispatch dance and are
                            # currently waiting on submissions to complete?
                            my @submissions = $state->submissions;
                            if (($state->dispatched || ($state->same_submissions_as && $state->same_submissions_as->dispatched)) && @submissions) {
                                warn "   had submissions\n" if $debug;
                                my @unfinished = VRPipe::Submission->get_column_values('id', { '_done' => 0, stepstate => $state->submission_search_id });
                                unless (@unfinished) {
                                    # check we're not the victim of a race condition and
                                    # that we still have submissions
                                    my $done = VRPipe::Submission->search({ '_done' => 1, stepstate => $state->submission_search_id });
                                    my $total = VRPipe::Submission->search({ stepstate => $state->submission_search_id });
                                    warn "   no unfinished subs; done $done vs total $total\n" if $debug;
                                    
                                    if ($total && $done == $total) {
                                        # run post_process
                                        my $pp_error;
                                        warn "   will post process\n" if $debug;
                                        $self->log_event("PipelineSetup->trigger, will post_process because all Submissions are done", dataelement => $element->id, stepstate => $state->id);
                                        eval { # (try catch not used because stupid perltidy is stupid)
                                            $pp_error = $step->post_process();
                                        };
                                        if ($@ && !$pp_error) {
                                            $pp_error = $@;
                                        }
                                        
                                        unless ($pp_error) {
                                            # we just completed all the submissions from a
                                            # previous parse
                                            $self->log_event("PipelineSetup->trigger found the post_process was successful", dataelement => $element->id, stepstate => $state->id);
                                            warn "   completed on pre-existing subs\n" if $debug;
                                            $self->_complete_state($step, $state, $step_number, $pipeline, $estate);
                                            $do_next = 1;
                                            
                                            # check to see if we've just completed the
                                            # whole setup
                                            if ($step_number == $num_steps) {
                                                if ($self->currently_complete) {
                                                    my $mt = VRPipe::MessageTracker->create(subject => "overall state of setup $setup_id");
                                                    my $num_states = $self->dataelementstates;
                                                    unless ($mt->already_sent("complete with $num_states elements")) {
                                                        $self->log_event("Completed setup with $num_states DataElements");
                                                        my $name = $self->name;
                                                        my $long = "\nTo remind yourself about this setup, do:\n\$ vrpipe-status --setup $setup_id\n\nTo get easy access to the output files, use vrpipe-output. eg:\n\$ vrpipe-output --setup $setup_id --output_with_input --basename_as_output\n\nIf this setup is now really complete (you won't be adding any more data to the datasource in future), please run:\n\$ vrpipe-setup --setup $setup_id --deactivate\n";
                                                        $self->_in_memory->log("Setup $setup_id ($name) has completed for $num_states DataElements", email_to => [$self->user], subject => "Setup $setup_id has completed", long_msg => $long);
                                                    }
                                                }
                                            }
                                            
                                            return;
                                        }
                                        else {
                                            # we warn instead of throw, because the step may
                                            # have discovered its output files are missing
                                            # and restarted itself
                                            $sm_error = "When trying to post process $error_ident after the submissions completed we hit the following error:\n$pp_error";
                                            warn "   $sm_error\n" if $debug;
                                            
                                            if ($pp_error =~ /start_over/) {
                                                warn "   - will redo since we just did a start_over\n" if $debug;
                                            }
                                            else {
                                                warn "   - this stepstate has permanent errors, will skip to the next elementstate\n" if $debug;
                                                $do_last = 1;
                                                return;
                                            }
                                        }
                                    }
                                    elsif (!$total) {
                                        $sm_error = "When dealing with what we thought were the completed submissions of $error_ident, the submissions vanished!";
                                        warn "   $sm_error; will die to restart transaction\n" if $debug;
                                        die $sm_error;
                                    }
                                    else {
                                        warn "   have subs but they're still running/failed - a handler should deal with them\n" if $debug;
                                    }
                                }
                                # else we have $unfinished unfinished submissions from a
                                # previous parse and are still running... unless we have
                                # the same submissions as something else, which might be
                                # for a setup that is no longer active...
                                elsif (!$error_message && $state->same_submissions_as) {
                                    my $other_state = $state->same_submissions_as;
                                    warn "   no error, and our state has the same subs as ", $other_state->id, "\n" if $debug;
                                    my $other_setup = $other_state->pipelinesetup;
                                    unless ($other_setup->active) {
                                        $sm_error = "The submissions for $error_ident were first created by setup " . $other_setup->id . ", but that setup is no longer active, so $setup_id is stalled!";
                                        warn "   - $sm_error\n" if $debug;
                                    }
                                    else {
                                        # we could also have same_submissions_as
                                        # something for a withdrawn dataelement, in
                                        # which case we'll also be stalled
                                        if ($other_state->dataelement->withdrawn) {
                                            $state->same_submissions_as(undef);
                                            $state->update;
                                            $self->log_event("PipelineSetup->trigger found that this StepState had the same submissions as for a withdrawn DataElement, so same_submissions_as was cleared", dataelement => $element->id, stepstate => $state->id);
                                            warn "   - same submissions as for a withdrawn dataelement, unset same_submissions_as\n" if $debug;
                                        }
                                        else {
                                            # else, is it safe to assume the submissions
                                            # of this other setup are really running?
                                            $self->debug("same submissions as");
                                            warn "   - assuming setup ", $other_setup->id, " will take care of things\n" if $debug;
                                        }
                                    }
                                }
                                else {
                                    warn "   had submissions and I guess they're running normally\n" if $debug;
                                    $do_last = 1;
                                    return;
                                }
                            }
                            else {
                                # this is the first time we're looking at this step for
                                # this data member for this pipelinesetup
                                my $parse_return;
                                warn "   no submissions, will parse\n" if $debug;
                                $self->log_event("PipelineSetup->trigger will parse", dataelement => $element->id, stepstate => $state->id);
                                try {
                                    $parse_return = $step->parse();
                                }
                                catch ($err) {
                                    $parse_return = $err;
                                }
                                
                                # if we're the primary for a block_and_skip job and have
                                # just been started over to get here, our stepoutputfile
                                # might have changed, so we need to update all the other
                                # stepstates that have same subs as us
                                unless ($parse_return) {
                                    warn "   - no parse error, will check output files are up-to-date on same_submissions_as states\n" if $debug;
                                    my %sof_file_to_key = map { $_->file->id => $_->output_key } $state->_output_files;
                                    foreach my $same (VRPipe::StepState->search({ same_submissions_as => $state->id })) {
                                        warn "    - fixing up state ", $same->id, "\n" if $debug;
                                        my %correct_sofs;
                                        foreach my $sof ($same->_output_files) {
                                            my $this_file = $sof->file->id;
                                            if (exists $sof_file_to_key{$this_file}) {
                                                unless ($sof->output_key eq $sof_file_to_key{$this_file}) {
                                                    $sof->output_key($sof_file_to_key{$this_file});
                                                    $sof->update;
                                                    $self->log_event("PipelineSetup->trigger after parse() for StepState " . $state->id . " found that the key for File $this_file changed to $sof_file_to_key{$this_file}, so StepOutputFile updated", dataelement => $same->dataelement->id, stepstate => $same->id);
                                                }
                                                $correct_sofs{$this_file} = 1;
                                            }
                                            else {
                                                $sof->delete;
                                                $self->log_event("PipelineSetup->trigger after parse() for StepState " . $state->id . " found that the $this_file was no longer an output, so StepOutputFile deleted", dataelement => $same->dataelement->id, stepstate => $same->id);
                                            }
                                        }
                                        
                                        while (my ($file, $key) = each %sof_file_to_key) {
                                            next if exists $correct_sofs{$file};
                                            VRPipe::StepOutputFile->create(stepstate => $same, file => $file, output_key => $key);
                                            $self->log_event("PipelineSetup->trigger after parse() for StepState " . $state->id . " found that $file was a new output, so StepOutputFile created", dataelement => $same->dataelement->id, stepstate => $same->id);
                                        }
                                        
                                        # also, set the dataelementstate to 0 and done to 0,
                                        # as if we had done a start_over on this state
                                        # (whilst avoiding the complications of really doing
                                        #  one)
                                        VRPipe::DataElementState->get(pipelinesetup => $same->pipelinesetup, dataelement => $same->dataelement, completed_steps => 0);
                                        $same->complete(0);
                                        $same->update;
                                        $self->log_event("PipelineSetup->trigger after parse() for StepState " . $state->id . " reset our complete to 0 since we had the same submissions as it", dataelement => $same->dataelement->id, stepstate => $same->id);
                                    }
                                }
                                
                                # $parse_return is 0 if we dispatched something, undef
                                # if we completed instantly and already ran post_process
                                # successfully, or is an error string
                                if ($parse_return) {
                                    $sm_error = "When trying to parse $error_ident we hit the following error:\n$parse_return";
                                    warn "   $sm_error\n" if $debug;
                                    
                                    #*** it would be better if the errors from
                                    # parse were classes so we could know exactly
                                    # the best thing to do here, but for now we'll
                                    # just parse the error message
                                    if ($parse_return =~ /start_over/) {
                                        warn "   - will redo since we just did a start_over\n" if $debug;
                                    }
                                    else {
                                        warn "   - this stepstate has permanent errors, will skip to the next elementstate\n" if $debug;
                                        $do_last = 1;
                                        return;
                                    }
                                }
                                elsif (!defined $parse_return) {
                                    $self->log_event("PipelineSetup->trigger called parse(), which dispatched nothing and completed instantly", dataelement => $element->id, stepstate => $state->id);
                                    $state->dispatched(1);
                                    $state->update;
                                    warn "   - instant complete after parse\n" if $debug;
                                    $self->_complete_state($step, $state, $step_number, $pipeline, $estate);
                                    $do_next = 1;
                                    return;
                                }
                                else {
                                    warn "   parse worked normally\n" if $debug;
                                    
                                    # if this step was supposed to have output
                                    # files, did it specify them?
                                    my $ohash = $state->output_files;
                                    my $odefs = $step->outputs_definition;
                                    foreach my $key (keys %$odefs) {
                                        next if exists $ohash->{$key};
                                        my $val = $odefs->{$key};
                                        next if $val->min_files == 0;
                                        $sm_error = "After parsing $error_ident we found that '$key' was defined as an output, yet no output file was made with that output_key";
                                        warn "   - $sm_error; will die to restart transaction\n" if $debug;
                                        die $sm_error;
                                    }
                                    
                                    my $dispatched = $step->dispatched();
                                    if (@$dispatched) {
                                        warn "   parse dispatched some jobs, will check if any other stepstates dispatched the same jobs\n" if $debug;
                                        
                                        # is there another stepstate that we already
                                        # made equivalent submissions for?
                                        my %other_states;
                                        my %this_cached_sameas;
                                        foreach my $arrayref (@$dispatched) {
                                            my ($cmd, undef, $job_args) = @$arrayref;
                                            
                                            if (defined $cached_sameas{$cmd}) {
                                                my ($ssid, $count) = @{ $cached_sameas{$cmd} };
                                                $other_states{$ssid} += $count;
                                                warn "   - stepstate $ssid did $count of our subs (cached)\n" if $debug;
                                                next;
                                            }
                                            
                                            my ($job) = VRPipe::Job->search({ dir => $output_root, $job_args ? (%{$job_args}) : (), cmd => $cmd });
                                            last unless $job;
                                            
                                            my @submissions = VRPipe::Submission->search({ job => $job->id, -or => [{ '_done' => 1 }, { 'dataelement.withdrawn' => 0 }] }, { prefetch => 'stepstate', join => { stepstate => 'dataelement' } });
                                            last unless @submissions;
                                            foreach my $sub (@submissions) {
                                                $other_states{ $sub->stepstate->id }++;
                                                $this_cached_sameas{ $sub->stepstate->id }->{$cmd}++;
                                                warn "   - stepstate ", $sub->stepstate->id, " did one of our subs (new)\n" if $debug;
                                            }
                                        }
                                        delete $other_states{ $state->id };
                                        
                                        my $same_as_us;
                                        my $needed_count = scalar @$dispatched;
                                        foreach my $other_id (sort { $a <=> $b } keys %other_states) {
                                            my $count = $other_states{$other_id};
                                            next unless $needed_count == $count;
                                            $same_as_us = $other_id;
                                            
                                            warn "   - stepstate $other_id had $count jobs the same as us\n" if $debug;
                                            
                                            while (my ($cmd, $count) = each %{ $this_cached_sameas{$other_id} }) {
                                                $cached_sameas{$cmd} = [$other_id, $count];
                                            }
                                            
                                            last;
                                        }
                                        
                                        if ($same_as_us) {
                                            # we just say that $state's
                                            # submissions are the same as the
                                            # other stepstate's
                                            $state->same_submissions_as($same_as_us);
                                            $state->update;
                                            $self->log_event("PipelineSetup->trigger called parse(), which dispatched the same submissions as StepState $same_as_us", dataelement => $element->id, stepstate => $state->id);
                                            
                                            warn "   - set same_submissions_as to $same_as_us, will redo and probably complete this step\n" if $debug;
                                            # (now we'll redo the loop; we probably
                                            # already completed this step)
                                        }
                                        else {
                                            warn "   - this is a new set of jobs, so will create jobs & submissions\n" if $debug;
                                            
                                            # in the case that this state has
                                            # been partially reset but all the
                                            # cmdlines have changed, we don't
                                            # want to create new subs for cmds
                                            # with output files that already
                                            # exist, and where the sub that made
                                            # those is done
                                            my %output_files_from_done_subs;
                                            foreach my $sub ($state->submissions()) {
                                                next unless $sub->_done;
                                                my $ofiles = $sub->job->output_files;
                                                if ($ofiles && @$ofiles) {
                                                    $output_files_from_done_subs{ join(',', sort map { $_->id } @$ofiles) } = 1;
                                                }
                                            }
                                            
                                            # create new submissions for the
                                            # relevant stepstate ($state may
                                            # have had start_over run on it,
                                            # which would have deleted its
                                            # same_submissions_as stepstate's
                                            # subs, and we want to create the
                                            # new submissions for that
                                            # same_submissions_as stepstate, not
                                            # for $state, hence the use of
                                            # $state->submission_search_id)
                                            $self->log_event("PipelineSetup->trigger called parse(), which dispatched new things", dataelement => $element->id, stepstate => $state->id);
                                            foreach my $arrayref (@$dispatched) {
                                                my ($cmd, $reqs, $job_args) = @$arrayref;
                                                
                                                if ($job_args && defined $job_args->{output_files}) {
                                                    if ($#{ $job_args->{output_files} } > 100) {
                                                        # protect us against
                                                        # job_args->output_files
                                                        # having too many values
                                                        # to fit in the db by
                                                        # just deleting it in
                                                        # that case: it's a
                                                        # nicety, not a
                                                        # necessity.
                                                        delete $job_args->{output_files};
                                                        undef $job_args unless keys %$job_args;
                                                    }
                                                    elsif (exists $output_files_from_done_subs{ join(',', sort map { $_->id } @{ $job_args->{output_files} }) }) {
                                                        # we've already made
                                                        # these files in an
                                                        # equivalent done sub
                                                        next;
                                                    }
                                                }
                                                
                                                # advertise that we're creating
                                                # subs with this reqs->id so
                                                # that if we create another one
                                                # soon vrpipe-server will put
                                                # them both in the same job
                                                # array
                                                $reqs->note('generating_subs', forget_after => 5);
                                                
                                                my $sub = VRPipe::Submission->create(job => VRPipe::Job->create(dir => $output_root, $job_args ? (%{$job_args}) : (), cmd => $cmd)->id, stepstate => $state->submission_search_id, requirements => $reqs->id);
                                                
                                                # because of globabl step limit
                                                # handling, we'll always need
                                                # vrpipe-server to look for all
                                                # subs and queue them to be
                                                # run, so there's not much value
                                                # in queuing $sub right now
                                                
                                                $self->log_event("PipelineSetup->trigger called parse(), and the dispatch created a new Submission", dataelement => $element->id, stepstate => $state->id, submission => $sub->id, job => $sub->job->id);
                                            }
                                            $state->dispatched(1);
                                            $state->update;
                                            warn "   - created new subs, will go to next elementstate\n" if $debug;
                                            $do_last = 1;
                                            return;
                                        }
                                    }
                                    else {
                                        # it is possible for a parse to result
                                        # in a different step being started over
                                        # because input files were missing
                                        warn "   - neither completed nor dispatched anything; it's possible a different step was started over due to missing input files?\n" if $debug;
                                    }
                                }
                            }
                        };
                        
                        eval { $self->do_transaction($transaction, "Failed to handle StepState " . $state->id . " during trigger()"); };
                        if ($@) {
                            $self->log_event($@, dataelement => $element->id, stepstate => $state->id);
                            warn "  main transaction failed: $@\n" if $debug;
                            $sm_error = $@;
                        }
                        else {
                            warn "  main transaction succeeded\n" if $debug;
                        }
                        
                        $do_next ||= 0;
                        $do_last ||= 0;
                        
                        if ($do_next) {
                            warn "  - going to next step member\n" if $debug;
                            next SSTATE;
                        }
                        elsif ($do_last) {
                            warn "  - last step member\n" if $debug;
                            $error_message ||= $sm_error if $sm_error;
                            last SSTATE;
                        }
                        
                        # if we got here we might have encountered some problem
                        # that fixed itself, so we'll try redoing the
                        # transaction
                        next;
                    }
                    
                    # if we got here we encountered problems more than
                    # $max_redos number of times; give up
                    $error_message ||= $sm_error;
                    my $warn_error = $sm_error || '[no error message]';
                    warn "  - giving up on this element state: $warn_error\n" if $debug;
                    last;
                }
                
                warn " finished with elementstate, will go on to next\n" if $debug;
                $estate->unlock unless $debug;
            }
        }
        
        my $warn_error = $error_message || '[no error message]';
        warn "trigger will return: $warn_error\n" if $debug;
        
        return $error_message;
    }
    
    method _complete_state (VRPipe::Step $step, VRPipe::StepState $state, Int $step_number, VRPipe::Pipeline $pipeline, VRPipe::DataElementState $estate) {
        my $transaction = sub {
            my $completed_steps = $estate->completed_steps;
            if ($step_number > $completed_steps) {
                $estate->completed_steps($step_number);
                $estate->update;
                $self->log_event("PipelineSetup->trigger found the DataElementState has now completed $step_number steps", dataelement => $estate->dataelement->id);
            }
            
            unless ($state->complete) {
                $self->log_event("PipelineSetup->trigger will now complete the StepState", dataelement => $estate->dataelement->id, stepstate => $state->id);
                
                # are there any behaviours to trigger?
                foreach my $behaviour (VRPipe::StepBehaviour->search({ pipeline => $pipeline->id, after_step => $step_number })) {
                    $self->log_event("PipelineSetup->trigger is triggering behaviour " . $behaviour->behaviour, dataelement => $estate->dataelement->id, stepstate => $state->id);
                    $behaviour->behave(data_element => $state->dataelement, pipeline_setup => $state->pipelinesetup);
                }
                
                unless ($state->same_submissions_as) {
                    # add to the StepStats
                    my %done_jobs;
                    foreach my $submission ($state->submissions) {
                        my $job = $submission->job;
                        next if $done_jobs{ $job->id };
                        VRPipe::StepStats->create(step => $step, pipelinesetup => $state->pipelinesetup, submission => $submission, memory => $job->peak_memory || 0, time => $job->wall_time);
                        $done_jobs{ $job->id } = 1;
                    }
                }
                
                $state->complete(1);
                $state->update;
            }
        };
        $self->do_transaction($transaction, "Failed to complete state for StepState " . $state->id);
    }
    
    method log_event (Str $msg, Bool :$record_stack?, Int :$dataelement = 0, Int :$stepstate = 0, Int :$submission = 0, Int :$job = 0) {
        #return unless $self->verbose; # at some point we'll only use this in testing?
        
        my $dt = DateTime->now();
        
        #*** PipelineSetupLogs don't always seem to get created, so first warn
        # what PipelineSetupLog->stringify would give us
        if ($self->verbose > 0) {
            my $friendly_dt = DateTime->from_epoch(epoch => $dt->epoch, time_zone => $local_timezone);
            my $date_str = "$friendly_dt";
            $date_str =~ s/T/ /;
            my $str = "$date_str | host $log_host | pid $$ | [ps " . $self->id;
            $str .= ", de $dataelement" if $dataelement;
            $str .= ", ss $stepstate"   if $stepstate;
            $str .= ", sub $submission" if $submission;
            $str .= ", job $job"        if $job;
            $str .= "] | " . $msg;
            chomp($str);
            warn $str, "\n";
            
            if ($record_stack) {
                warn $self->stack_trace, "\n";
            }
        }
        
        # the msg can't be longer than 65535 chars or it will fail to be stored
        # in the db. And actually, if its longer than a 1000 it's crazy, so we
        # truncate
        if (length($msg) > 1000) {
            $msg = substr($msg, 0, 1000) . ' ...[truncated]...';
        }
        
        return VRPipe::PipelineSetupLog->create(
            ps_id   => $self->id,
            message => $msg,
            date    => $dt,      # dates in the database are not in local time
            host    => $log_host,
            pid     => $$,
            $record_stack ? (stack => $self->stack_trace) : (),
            de_id        => $dataelement,
            ss_id        => $stepstate,
            sub_id       => $submission,
            job_id       => $job,
            blind_create => 1            # without this we hit deadlocks constantly, and we always want to create a new row anyway
        );
    }
    
    method logs (Str :$like?, Int :$dataelement?, Int :$stepstate?, Int :$submission?, Int :$job?, Int :$pid?, Str :$host?, Bool :$include_undefined?, Bool :$paged?, Str :$from?, Str :$to?) {
        my @date_args = ();
        if ($from || $to) {
            my $parser = DateTime::Format::Natural->new(time_zone => $local_timezone);
            my @dts;
            foreach my $dstr ($from, $to) {
                unless ($dstr) {
                    push(@dts, undef);
                    next;
                }
                
                # allow copy/pastes of stringified DateTimes
                $dstr =~ s/(\d)T(\d)/$1 $2/;
                
                my $dt = $parser->parse_datetime($dstr);
                if ($parser->success) {
                    push(@dts, DateTime->from_epoch(epoch => $dt->epoch));
                }
                else {
                    warn $parser->error, "\n";
                    push(@dts, undef);
                }
            }
            
            if ($dts[0] && !$dts[1]) {
                push(@date_args, date => { '>=' => $dts[0] });
            }
            elsif (!$dts[0] && $dts[1]) {
                push(@date_args, date => { '<=' => $dts[1] });
            }
            elsif ($dts[0] && $dts[1]) {
                push(@date_args, -and => [{ date => { '>=' => $dts[0] } }, { date => { '<=' => $dts[1] } }]);
            }
        }
        
        my $method = $paged ? 'search_paged' : 'search';
        return VRPipe::PipelineSetupLog->$method({
                ps_id => $self->id,
                $like ? (message => { like => $like }) : (),
                $dataelement ? (de_id  => $include_undefined ? { -in => [0, $dataelement] } : $dataelement) : (),
                $stepstate   ? (ss_id  => $include_undefined ? { -in => [0, $stepstate] }   : $stepstate)   : (),
                $submission  ? (sub_id => $include_undefined ? { -in => [0, $submission] }  : $submission)  : (),
                $job         ? (job_id => $include_undefined ? { -in => [0, $job] }         : $job)         : (),
                $pid  ? (pid  => $pid)  : (),
                $host ? (host => $host) : (),
                @date_args
            }
        );
    }
}

1;
