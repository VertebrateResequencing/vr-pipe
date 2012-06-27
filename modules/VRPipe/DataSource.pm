=head1 NAME

VRPipe::DataSource - define pipeline input from arbitrary sources

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A DataSource defines the inputs you want to supply to one or more
L<VRPipe::PipelineSetup>s. There are a number of different types of DataSource.
Each type understands a different kind of source. For example the 'fofn' type
understands text files that list file paths. The 'vrtrack' type understands
the contents of a VRTrack database.

Each type has one or more methods, and each method can have options. For
example, you could use the 'fofn' type with the 'all' method and no options,
and give it a 'source' of a file of filenames. That DataSource would then
define that each file in your fofn should be an independent input to any
PipelineSetup using this DataSource.

For example, consider a L<VRPipe::Pipeline> with 2 steps. The first step takes
a txt file as input, copies it, and appends a line of text to the copy, then
outputs the new text file. The second step takes the output of the first, makes
a copy of it, and removes the first line before outputting this altered version.

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
is input into the first step. In parallel the code of the first step us run
on the inputs, and in the first working directory a file 11 lines long is made,
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
    use VRPipe::DataSourceFactory;
    
    has 'type' => (is => 'rw',
                   isa => Varchar[64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'method' => (is => 'rw',
                     isa => Varchar[64],
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'source' => (is => 'rw',
                     isa => Text,
                     coerce => 1,
                     traits => ['VRPipe::Persistent::Attributes'],
                     is_key => 1);
    
    has 'options' => (is => 'rw',
                      isa => 'HashRef',
                      traits => ['VRPipe::Persistent::Attributes'],
                      default => sub { {} },
                      allow_key_to_default => 1,
                      is_key => 1);
    
    has '_changed_marker' => (is => 'rw',
                             isa => Varchar[255],
                             traits => ['VRPipe::Persistent::Attributes'],
                             is_nullable => 1);
    
    has '_source_instance' => (is => 'rw',
                               isa => 'Defined',
                               lazy => 1,
                               builder => '_build_source',
                               handles => [qw(description source_description method_description)]);
    
    method _build_source {
        my $changed_marker = $self->_changed_marker;
        return VRPipe::DataSourceFactory->create($self->type, {method => $self->method,
                                                               source => $self->source,
                                                               options => $self->options,
                                                               $changed_marker ? ('_changed_marker' => $changed_marker) : (),
                                                               '_datasource_id' => $self->id});
    }
    
    __PACKAGE__->make_persistent(has_many => [elements => 'VRPipe::DataElement']);
    
    around elements {
        # we don't return $self->$orig because that returns all associated
        # elements; we need to first create all elements, and then only return
        # elements that are still "current" (haven't been deleted in the source)
        $self->_prepare_elements_and_states || return;
        
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('DataElement')->search({ datasource => $self->id, withdrawn => 0 });
        
        my @elements;
        while (my $element = $rs->next) {
            push(@elements, $element);
        }
        
        return \@elements;
    }
    
    method incomplete_element_states (VRPipe::PipelineSetup $setup, Int $limit?) {
        $self->_prepare_elements_and_states || return;
        
        my $pipeline = $setup->pipeline;
        my $num_steps = $pipeline->step_members;
        
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('DataElementState')->search({ pipelinesetup => $setup->id, completed_steps => {'<', $num_steps}, 'dataelement.withdrawn' => 0 },
                                                                { join => 'dataelement', $limit ? (rows => $limit) : () });
        
        my @incomplete;
        while (my $state = $rs->next) {
            push(@incomplete, $state);
        }
        
        return \@incomplete;
    }
    
    method _prepare_elements_and_states {
        my $source = $self->_source_instance || return;
        
        my $schema = $self->result_source->schema;
        my $rs = $schema->resultset('PipelineSetup')->search({ datasource => $self->id });
        my @setups;
        while (my $ps = $rs->next) {
            push(@setups, $ps);
        }
        
        if ($source->_has_changed) {
            my $elements = $source->_get_elements;
            
            my %current_elements;
            foreach my $element (@$elements) {
                $current_elements{$element->id} = 1;
                foreach my $setup (@setups) {
                    VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $element);
                }
            }
            
            # withdrawn any elements that are no longer in the datasource
            my $schema = $self->result_source->schema;
            my $rs = $schema->resultset('DataElement')->search({ datasource => $self->id, withdrawn => 0 });
            while (my $element = $rs->next) {
                unless (exists $current_elements{$element->id}) {
                    $element->withdrawn(1);
                    $element->update;
                }
            }
            
            $self->_changed_marker($source->_changed_marker);
            $self->update;
        }
        else {
            # check that a new pipelinesetup hasn't been created since the source
            # last changed
            my $expected_count = $schema->resultset('DataElement')->count({ datasource => $self->id });
            foreach my $setup (@setups) {
                my $count = $schema->resultset('DataElementState')->count({ pipelinesetup => $setup->id });
                if ($count < $expected_count) {
                    $rs = $schema->resultset('DataElement')->search({ datasource => $self->id });
                    my @elements;
                    while (my $element = $rs->next) {
                        VRPipe::DataElementState->get(pipelinesetup => $setup, dataelement => $element);
                    }
                }
            }
        }
        
        return 1;
    }
}

1;
