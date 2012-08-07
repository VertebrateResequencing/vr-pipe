
=head1 NAME

VRPipe::DataSource - define pipeline input from arbitrary sources

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A DataSource defines the inputs you want to supply to one or more
L<VRPipe::PipelineSetup>s. There are a number of different types of DataSource.
Each type understands a different kind of source. For example the 'fofn' type
understands text files that list file paths. The 'vrtrack' type understands the
contents of a VRTrack database.

Each type has one or more methods, and each method can have options. For
example, you could use the 'fofn' type with the 'all' method and no options,
and give it a 'source' of a file of filenames. That DataSource would then
define that each file in your fofn should be an independent input to any
PipelineSetup using this DataSource.

For example, consider a L<VRPipe::Pipeline> with 2 steps. The first step takes
a txt file as input, copies it, and appends a line of text to the copy, then
outputs the new text file. The second step takes the output of the first, makes
a copy of it, and removes the first line before outputting this altered
version.

Now let's say we make a DataSource of type 'fofn', method 'all', source
'/abs/that/to/my.fofn', where the contents of /abs/that/to/my.fofn are:
    
    /abs/path/to/file.1
    /abs/path/to/file.2

/abs/path/to/file.1 is a text file 10 lines long, and /abs/path/to/file.2 is a
text file 5 lines long.

This DataSource will be comprised of 2 L<VRPipe::DataElement>s, one for
/abs/path/to/file.1, and the other for /abs/path/to/file.2

Finally, we make a PipelineSetup that uses the previously mentioned Pipeline
and DataSource. What will happen is this:

In a certain working directory, /abs/path/to/file.1 is input into the first
step. At the same time, in a different working directory, /abs/path/to/file.2
is input into the first step. In parallel the code of the first step us run on
the inputs, and in the first working directory a file 11 lines long is made,
and in the second working directory a file of 6 lines long is made.

Soon after the 11 line long file is made, it is passed to the second step.
Independently (but probably around the same time in this example), soon after
the 6 line long file is made, it is passed to the second step. In a third
working directory the 11 line long file is processed to generate a 10 line long
file, and (probably in parallel) in a fourth working directory the 6 line long
file is processed to generate a 5 line long file.

At this point the pipeline is complete, and you can find the 2 output files in
the third and fourth working directories and see which input file they
correspond to by using the B<vrpipe-output> script.

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

class VRPipe::DataSource extends VRPipe::Persistent {
    use DateTime;
    use VRPipe::DataSourceFactory;
    
    has 'type' => (is     => 'rw',
                   isa    => Varchar [64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'method' => (is     => 'rw',
                     isa    => Varchar [64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'source' => (is     => 'rw',
                     isa    => Text,
                     coerce => 1,
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'options' => (is                   => 'rw',
                      isa                  => 'HashRef',
                      traits               => ['VRPipe::Persistent::Attributes'],
                      default              => sub { {} },
                      allow_key_to_default => 1,
                      is_key               => 1);
    
    has '_lock' => (is          => 'rw',
                    isa         => Datetime,
                    traits      => ['VRPipe::Persistent::Attributes'],
                    is_nullable => 1);
    
    has '_changed_marker' => (is          => 'rw',
                              isa         => Varchar [255],
                              traits      => ['VRPipe::Persistent::Attributes'],
                              is_nullable => 1);
    
    has '_source_instance' => (is      => 'rw',
                               isa     => 'Defined',
                               lazy    => 1,
                               builder => '_build_source',
                               handles => [qw(description source_description method_description)]);
    
    method _build_source {
        my $changed_marker = $self->_changed_marker;
        return
          VRPipe::DataSourceFactory->create($self->type,
                                            {  method  => $self->method,
                                               source  => $self->source,
                                               options => $self->options,
                                               $changed_marker ? ('_changed_marker' => $changed_marker) : (),
                                               '_datasource_id' => $self->id });
    }
    
    __PACKAGE__->make_persistent(has_many => [elements => 'VRPipe::DataElement']);
    
    around elements {
        my (undef, undef, $status_messages) = @_;
        # we don't return $self->$orig because that returns all associated
        # elements; we need to first create all elements, and then only return
        # elements that are still "current" (haven't been deleted in the source)
        $self->_prepare_elements_and_states($status_messages) || return;
        return VRPipe::DataElement->search_paged({ datasource => $self->id, withdrawn => 0 });
    }
    
    method incomplete_element_states (VRPipe::PipelineSetup $setup, Bool :$include_withdrawn = 0, Bool :$prepare = 1) {
        if ($prepare) {
            $self->_prepare_elements_and_states || return;
        }
        
        my $pipeline  = $setup->pipeline;
        my $num_steps = $pipeline->step_members;
        
        return VRPipe::DataElementState->search_paged({ pipelinesetup => $setup->id, completed_steps => { '<', $num_steps }, $include_withdrawn ? () : ('dataelement.withdrawn' => 0) }, { prefetch => 'dataelement' });
    }
    
    method _prepare_elements_and_states (Bool $status_messages = 0) {
        my $source = $self->_source_instance || return;
        
        my @setup_ids = VRPipe::PipelineSetup->get_column_values('id', { datasource => $self->id });
        
        if ($source->_has_changed) {
            # we must not go through and update the dataelements more than
            # once simultaneously, and we must not return elements in a
            # partially updated state, so we lock/block at this point
            my $block    = 0;
            my $continue = 1;
            do {
                $self->reselect_values_from_db;
                my $transaction = sub {
                    my $lock_time = $self->_lock;
                    # check that the process that got the lock is still running,
                    # otherwise ignore the lock
                    if ($lock_time) {
                        my $elapsed = time() - $lock_time->epoch;
                        if ($elapsed > 60) {
                            undef $lock_time;
                        }
                    }
                    
                    # if some other process has a recent lock, we will block
                    if ($lock_time) {
                        $block = 1;
                        warn "another process is updating dataelements, please wait...\n" if $status_messages;
                        sleep(2);
                        return;
                    }
                    
                    # if we had been blocking and now there is no more lock,
                    # likely that the datasource is now up to date and we don't
                    # have to do anything
                    if ($block) {
                        $block    = 0;
                        $continue = $source->_has_changed;
                    }
                    
                    if ($continue) {
                        # get the lock for ourselves
                        $self->_lock(DateTime->now());
                        $self->update;
                    }
                };
                $self->do_transaction($transaction, "DataSource lock/block failed");
            } while ($block);
            return unless $continue;
            
            # we have a lock, but can't risk our process getting killed and the
            # lock being left open, so we have to re-claim the lock every 15s
            # so that the above blocking code doesn't ignore the lock. We do
            # this by spawning a lock process
            my $my_pid   = $$;
            my $lock_pid = fork();
            if (!defined $lock_pid) {
                $self->throw("attempt to fork for lock failed: $!");
            }
            elsif ($lock_pid == 0) {
                # child, initiate a lock that will end when the parent stops
                # running
                sleep(2);
                while (1) {
                    kill(0, $my_pid) || last;
                    warn "updating dataelements, please wait...\n" if $status_messages;
                    $self->_lock(DateTime->now());
                    $self->update;
                    $self->disconnect;
                    sleep 15;
                }
                exit(0);
            }
            
            # we need to create a DataElementState for each setup using this
            # datasource for each new dataelement. To minimise time wasted
            # trying to create DES when they already exist, before getting new
            # elements we find out what the most recent element id is:
            my ($most_recent_element_id) = VRPipe::DataElement->get_column_values('id', { datasource => $self->id }, { order_by => { -desc => ['id'] }, rows => 1 });
            
            # now get the source to create new dataelements:
            $source->_generate_elements;
            
            # now page through dataelements with a higher id than previous most
            # recent:
            my @des_args;
            my $pager = VRPipe::DataElement->get_column_values_paged('id', { datasource => $self->id, $most_recent_element_id ? (id => { '>' => $most_recent_element_id }) : () });
            while (my $eids = $pager->next) {
                foreach my $eid (@$eids) {
                    foreach my $setup_id (@setup_ids) {
                        push(@des_args, { pipelinesetup => $setup_id, dataelement => $eid });
                    }
                }
            }
            VRPipe::DataElementState->bulk_create_or_update(@des_args) if @des_args;
            
            # we're done, so update changed marker
            $self->_changed_marker($source->_changed_marker);
            $self->update;
            
            # cleanup the lock process and set _lock() to 1 minute ago (we can't
            # undef it, but doing this gives the same result)
            kill(9, $lock_pid);
            waitpid($lock_pid, 0);
            $self->_lock(DateTime->from_epoch(epoch => time() - 60));
            $self->update;
        }
        else {
            # check that a new pipelinesetup hasn't been created since the source
            # last changed
            my $expected_count = VRPipe::DataElement->search({ datasource => $self->id });
            
            my @des_args;
            foreach my $setup_id (@setup_ids) {
                my $count = VRPipe::DataElementState->search({ pipelinesetup => $setup_id });
                
                if ($count < $expected_count) {
                    my $pager = VRPipe::DataElement->get_column_values_paged('id', { datasource => $self->id });
                    while (my $dataelement_ids = $pager->next) {
                        foreach my $eid (@$dataelement_ids) {
                            push(@des_args, { pipelinesetup => $setup_id, dataelement => $eid });
                        }
                    }
                }
            }
            VRPipe::DataElementState->bulk_create_or_update(@des_args) if @des_args;
        }
        
        return 1;
    }
}

1;
