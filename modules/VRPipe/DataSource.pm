
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
    
    has 'type' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'method' => (
        is     => 'rw',
        isa    => Varchar [64],
        traits => ['VRPipe::Persistent::Attributes'],
        is_key => 1
    );
    
    has 'source' => (
        is     => 'rw',
        isa    => Text,
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
    
    has '_changed_marker' => (
        is          => 'rw',
        isa         => Varchar [255],
        traits      => ['VRPipe::Persistent::Attributes'],
        is_nullable => 1
    );
    
    has '_source_instance' => (
        is      => 'rw',
        isa     => 'Defined',
        lazy    => 1,
        builder => '_build_source',
        handles => [qw(description source_description method_description)]
    );
    
    has 'debug' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    method _build_source {
        my $changed_marker = $self->_changed_marker;
        return VRPipe::DataSourceFactory->create(
            $self->type,
            {
                method  => $self->method,
                source  => $self->source,
                options => $self->options,
                $changed_marker ? ('_changed_marker' => $changed_marker) : (),
                '_datasource_id' => $self->id,
                debug            => $self->debug
            }
        );
    }
    
    __PACKAGE__->make_persistent(has_many => [elements => 'VRPipe::DataElement']);
    
    around elements {
        my (undef, undef, $status_messages, $debug) = @_;
        # we don't return $self->$orig because that returns all associated
        # elements; we need to first create all elements, and then only return
        # elements that are still "current" (haven't been deleted in the source)
        $self->_prepare_elements_and_states($status_messages, $debug || $self->debug);
        return VRPipe::DataElement->search_paged({ datasource => $self->id, withdrawn => 0 });
    }
    
    method incomplete_element_states (VRPipe::PipelineSetup $setup, Bool :$include_withdrawn = 0, Bool :$prepare = 1, Bool :$only_not_started = 0, Bool :$debug = 0) {
        if ($prepare) {
            $self->_prepare_elements_and_states(0, $debug || $self->debug);
        }
        
        my $pipeline  = $setup->pipeline;
        my $num_steps = $pipeline->step_members;
        
        return VRPipe::DataElementState->search_paged({ pipelinesetup => $setup->id, completed_steps => $only_not_started ? 0 : { '<', $num_steps }, $include_withdrawn ? () : ('dataelement.withdrawn' => 0) }, { prefetch => 'dataelement' });
    }
    
    # for critical columns that get updated by multiple processes, ensure we
    # always get the latest value from db
    around _changed_marker (Str $marker?) {
        if ($marker) {
            return $self->$orig($marker);
        }
        $self->reselect_values_from_db;
        return $self->$orig;
    }
    
    method _prepare_elements_and_states (Bool $status_messages?, Bool $debug?) {
        my $source = $self->_source_instance || return;
        $debug ||= $self->debug;
        unless (defined $status_messages) {
            $status_messages = $self->verbose > 0 ? 1 : 0;
        }
        
        # we must not go through and update the dataelements more than once
        # simultaneously, and we must not return elements in a partially updated
        # state, so we lock/block at this point
        warn "will check the datasource, first blocking until locked\n" if $debug;
        my $blocked = $self->block_until_locked(debug => $debug);
        if (1 || $blocked) { #*** I think it's an artifact of the test in DataSource.t that I have to force this reselect
            $self->reselect_values_from_db;
            my $current_marker = $self->_changed_marker;
            if ($current_marker) {
                $source->_changed_marker($current_marker);
            }
        }
        warn "locked the datasource\n" if $debug;
        
        my @setup_ids = VRPipe::PipelineSetup->get_column_values('id', { datasource => $self->id });
        $source->debug($debug);
        
        my $generated_elements = 0;
        if ($source->_has_changed) {
            warn "the datasource source has changed\n" if $debug;
            # we need to create a DataElementState for each setup using this
            # datasource for each new dataelement. To minimise time wasted
            # trying to create DES when they already exist, before getting new
            # elements we find out what the most recent element id is:
            my ($most_recent_element_id) = VRPipe::DataElement->get_column_values('id', { datasource => $self->id }, { order_by => { -desc => ['id'] }, rows => 1 });
            
            # now get the source to create new dataelements:
            $source->_generate_elements;
            $generated_elements = 1;
            
            # now page through dataelements with a higher id than previous most
            # recent:
            my (@des_args, @enqueue_args);
            my $pager = VRPipe::DataElement->get_column_values_paged('id', { datasource => $self->id, $most_recent_element_id ? (id => { '>' => $most_recent_element_id }) : () });
            while (my $eids = $pager->next) {
                foreach my $eid (@$eids) {
                    foreach my $setup_id (@setup_ids) {
                        push(@des_args, { pipelinesetup => $setup_id, dataelement => $eid });
                        
                        # bulk_create_or_update call below alters @des_args, so
                        # we push to another array for queueing purposes
                        push(@enqueue_args, [$setup_id, $eid]);
                    }
                }
            }
            warn "generated ", scalar(@des_args), " new dataelements, will now create the states for them\n" if $debug;
            VRPipe::DataElementState->bulk_create_or_update(@des_args) if @des_args;
            
            # queue up these new de to be triggered (vrpipe-server will do
            # something with these); we don't queue des ids because we don't
            # know them, and that's because it's much faster for the above call
            # to be in void context
            my $im = $self->_in_memory;
            foreach my $args (@enqueue_args) {
                $im->enqueue('trigger', "$args->[0]:$args->[1]");
            }
            
            # we're done, so update changed marker
            my $cm = $source->_changed_marker;
            $self->_changed_marker($cm);
            $self->update;
            warn "set _changed_marker to $cm\n" if $debug;
        }
        else {
            warn "the datasource source hasn't changed\n" if $debug;
            # check that a new pipelinesetup hasn't been created since the source
            # last changed
            my $expected_count = VRPipe::DataElement->search({ datasource => $self->id });
            
            my (@des_args, @enqueue_args);
            foreach my $setup_id (@setup_ids) {
                my $count = VRPipe::DataElementState->search({ pipelinesetup => $setup_id });
                
                if ($count < $expected_count) {
                    my $pager = VRPipe::DataElement->get_column_values_paged('id', { datasource => $self->id });
                    while (my $dataelement_ids = $pager->next) {
                        foreach my $eid (@$dataelement_ids) {
                            push(@des_args, { pipelinesetup => $setup_id, dataelement => $eid });
                            push(@enqueue_args, [$setup_id, $eid]);
                        }
                    }
                }
            }
            
            if (@des_args) {
                warn "new setups were found, will create ", scalar(@des_args), " dataelementstates for them\n" if $debug;
                VRPipe::DataElementState->bulk_create_or_update(@des_args);
                
                my $im = $self->_in_memory;
                foreach my $args (@enqueue_args) {
                    $im->enqueue('trigger', "$args->[0]:$args->[1]");
                }
            }
        }
        $self->unlock;
        
        warn "all done, will return generated dataelements\n" if $debug;
        return $generated_elements;
    }
}

1;
