
=head1 NAME

VRPipe::Persistent::Pager - easily work with paged search results

=head1 SYNOPSIS
    
    use VRPipe::Persistent::Schema;
    
    my $pager = VRPipe::StepStats->search_paged({...});
    
    while (my $stepstats = $pager->next) {
	# $stepstats is an array ref of VRPipe::StepStats instances
    }
    
    # create your own Pager with your own ResultSet
    my $rs = VRPipe::StepStats->search_rs({...}); # you don't need row/page attributes
    $pager = VRPipe::Persistent::Pager->new(resultset => $rs, rows_per_page => 1000);

=head1 DESCRIPTION

When you want to work with many rows from the database you might use up too
much memory to get and store them all in one go with search(), but getting a
ResultSet with search_rs() and doing an $rs->next loop would be too slow since
that selects out each row one at a time with each next call.

Instead, methods like search_paged() create ResultSets that are paged and
create an instance of this class with that ResultSet. This class then makes it
easy to access a page's worth of the search results at a time.

If memory is sufficient, you should try 10000 rows per page for good
efficiency. The default is only 1000, which could be less than half the speed
of 10000.

NB: if what you are searching for changes between calls of next(), the query
will be restarted and you will be on page 1 again (possibly going through rows
you already dealt with, depending on what changed and what your query was). To
help avoid near-infinite loops if paging whilst new rows are added to the table
you're searching on, an extra search term is automatically applied that limits
the id to the last id when search_paged as first called.

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

class VRPipe::Persistent::Pager {
    has resultset => (
        is       => 'ro',
        isa      => 'DBIx::Class::ResultSet',
        required => 1,
        writer   => '_modify_resultset'
    );
    
    has rows_per_page => (
        is      => 'ro',
        default => 1000
    );
    
    has result_method => (
        is      => 'ro',
        isa     => 'Str|CodeRef',
        default => 'all'
    );
    
    has _pager => (
        is      => 'rw',
        isa     => 'Data::Page',
        lazy    => 1,
        builder => '_build_pager',
        handles => [qw(current_page last_page total_entries)]
    );
    
    has _pages_done => (
        is      => 'rw',
        default => 0
    );
    
    method _build_pager {
        my $rs = $self->resultset;
        
        # because another process could create new table rows which would change
        # our results every next(), we avoid unecessary (seemingly
        # psuedo-infinite) restarts by first getting the most recently created
        # matching row and then searching for ids less than that.
        my @id_limit = ();
        unless (defined $rs->{cond}->{'me.id'}) {
            my @refs = $rs->search_rs({}, { columns => ['me.id'] })->cursor->all; # avoid 'order_by => { -desc => 'me.id' }, rows => 1,' because it may be faster to send more data and sort in perl than it is to sort in the db
            if (@refs && $refs[0]->[0]) {
                my $max_id = 0;
                foreach my $ref (@refs) {
                    if ($ref->[0] > $max_id) {
                        $max_id = $ref->[0];
                    }
                }
                @id_limit = ('me.id' => { '<=' => $max_id });
            }
        }
        
        $rs = $rs->search({@id_limit}, { rows => $self->rows_per_page, page => 1 });
        $self->_modify_resultset($rs);
        return $rs->pager;
    }
    
    method next (Bool :$no_resetting?) {
        my $current_entries = $self->total_entries; # we must build our pager first, which alters resultset()
        my $rs              = $self->resultset;
        unless ($no_resetting) {
            my $new_rs = $rs->search({}, { rows => $self->rows_per_page, page => 1 });
            if ($new_rs->pager->total_entries != $current_entries) {
                $rs = $new_rs;
                $self->_modify_resultset($rs);
                $self->_pager($rs->pager);
                $self->_pages_done(0);
            }
        }
        
        my $last_page = $self->last_page;
        my $next_page = $self->_pages_done + 1;
        return if $next_page > $last_page;
        $rs = $rs->page($next_page);
        
        my $result_method = $self->result_method;
        my $rows;
        if (ref($result_method)) {
            my @results = &$result_method($rs);
            if (@results == 1 && ref($results[0]) eq 'ARRAY') {
                $rows = $results[0];
            }
            else {
                $rows = \@results;
            }
        }
        else {
            $rows = [$rs->$result_method];
        }
        
        $self->_pages_done($next_page);
        
        return $rows;
    }
    
    method reset {
        $self->_pages_done(0);
        $self->current_page(1);
    }
}

1;
