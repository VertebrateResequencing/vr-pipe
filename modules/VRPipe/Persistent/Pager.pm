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

Instead, methods like search_paged() create ResultSets that are paged and create
an instance of this class with that ResultSet. This class then makes it easy
to access a page's worth of the search results at a time.

If memory is sufficient, you should try 10000 rows per page for good efficiency.
The default is only 1000, which could be less than half the speed of 10000.

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
    has resultset => (is => 'ro',
                      isa => 'DBIx::Class::ResultSet',
                      required => 1,
		      writer => '_modify_resultset');
    
    has rows_per_page => (is => 'ro',
			  default => 1000);
    
    has _pager => (is => 'ro',
		   isa => 'Data::Page',
		   lazy => 1,
		   builder => '_build_pager',
		   handles => [qw(current_page last_page)]);
    
    has _pages_done => (is => 'rw',
			default => 0);
    
    method _build_pager {
	my $rs = $self->resultset;
	$rs = $rs->search({}, { rows => $self->rows_per_page, page => 1 });
	$self->_modify_resultset($rs);
	return $rs->pager;
    }
    
    method next {
	my $last_page = $self->last_page;
	my $next_page = $self->_pages_done + 1;
	return if $next_page > $last_page;
	
	my $rows = [$self->resultset->page($next_page)->all];
	
	$self->_pages_done($next_page);
	
	return $rows;
    }
    
    method reset {
	$self->_pages_done(0);
	$self->current_page(1);
    }
}

1;
