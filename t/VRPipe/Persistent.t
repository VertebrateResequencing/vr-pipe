#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 30;
    
    use_ok('VRPipe::Persistent');
    use_ok('t::VRPipe::Persistent::Schema');
    
    use Test::DBIx::Class {
        schema_class => 't::VRPipe::Persistent::Schema',
        force_drop_table => 1,
        keep_db => 1,
        connect_info => ["dbi:mysql:host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT};database=vertres_dbixclass;", $ENV{VRTRACK_RW_USER}, $ENV{VRTRACK_PASSWORD}] #['dbi:SQLite:dbname=:memory:','','']
    }; #, 't::VRPipe::Artist' => {-as => 'Artist'}, 't::VRPipe::CD' => {-as => 'CD'}, 't::VRPipe::Track' => {-as => 'Track'}
}

# some quick basic tests
ok my $resultset = ResultSet('Artist'), "Yeah, some Artists in the database";

ok my $artist = $resultset->create({name => 'Barton'}), 'Created Barton';
ok $artist = $resultset->create({name => 'Footon'}), 'Created Footon';

ok my $footon = $resultset->search({'name' => 'Footon'})->first, 'Searched for Footon';

is_fields [qw/artistid name age/], $footon, [2, 'Footon', 99], 'Footon has the expected fields';


# update behaviour
$artist->name('Gerton');
is $artist->name, 'Gerton', 'able to change name to Gerton';
is_deeply [$artist->get_dirty_columns], [name => 'Gerton'], 'name column marked as dirty after change, prior to update';
$artist->update;
is_deeply [$artist->get_dirty_columns], [], 'no more dirty columns after an update';
undef $artist;
ok $artist = $resultset->search({'name' => 'Gerton'})->first, 'Got Gerton back';
is $artist->name, 'Gerton', 'the name is Gerton';
is_deeply [$artist->id], [2], 'Gerton has the correct id';
is $artist->name('Larton'), 'Larton', 'able to change name again to Larton';
undef $artist;
ok $artist = $resultset->search({'artistid' => 2})->first, 'Got artist 2 back after undef';
is $artist->name, 'Larton', 'the name is Larton even after object destruction with no explicit update';


# constraint checking
throws_ok { $artist->name('12345678901234567890123456789012345678901234567890123456789012345') } qr/too long/, 'giving a name that is too long causes a throw';
is $artist->name, 'Larton', 'and the name is unchanged';


# throw behaviour
throws_ok { $resultset->search({'named' => 'Footon'})->first, 'Searched for Footon'; } 'DBIx::Class::Exception', 'throws given invalid column name to search on';


# try out our special get method
ok my $bob = t::VRPipe::Artist->get(schema => Schema, name => 'Bob'), 'created Bob using get()';
is_fields [qw/artistid name age/], $bob, [3, 'Bob', 99], 'Bob has the expected fields';
undef $bob;
ok $bob = t::VRPipe::Artist->get(schema => Schema, name => 'Bob'), 'got Bob using get()';
is_fields [qw/artistid name age/], $bob, [3, 'Bob', 99], 'Bob still has the expected fields';


# relationship testing from dbic example code
reset_schema; # this counts as a test
{
    # populate the database
    my @artists = (['Michael Jackson'], ['Eminem']);
    Schema->populate('Artist', [
        [qw/name/],
        @artists,
    ]);
    
    my %albums = (
        'Thriller' => 'Michael Jackson',
        'Bad' => 'Michael Jackson',
        'The Marshall Mathers LP' => 'Eminem',
    );
    
    my @cds;
    foreach my $lp (keys %albums) {
        my $artist = Schema->resultset('Artist')->find({
            name => $albums{$lp}
        });
        push @cds, [$lp, $artist->id];
    }
    
    Schema->populate('CD', [
        [qw/title artist/],
        @cds,
    ]);
    
    my %tracks = (
        'Beat It'         => 'Thriller',
        'Billie Jean'     => 'Thriller',
        'Dirty Diana'     => 'Bad',
        'Smooth Criminal' => 'Bad',
        'Leave Me Alone'  => 'Bad',
        'Stan'            => 'The Marshall Mathers LP',
        'The Way I Am'    => 'The Marshall Mathers LP',
    );
    
    my @tracks;
    foreach my $track (keys %tracks) {
        my $cdname = Schema->resultset('CD')->find({
            title => $tracks{$track},
        });
        push @tracks, [$cdname->id, $track];
    }
    
    Schema->populate('Track',[
        [qw/cd title/],
        @tracks,
    ]);
    
    # did we create the correct tables in the db?
    #... difficult to test given configurable dbm
    
    # run tests
    is_deeply [get_tracks_by_cd('Bad')], ['Leave Me Alone', 'Smooth Criminal', 'Dirty Diana'], 'could get tracks by cd';
    is_deeply [get_tracks_by_artist('Michael Jackson')], ['Billie Jean', 'Beat It', 'Leave Me Alone', 'Smooth Criminal', 'Dirty Diana'], 'could get tracks by artist';
    
    is get_cd_by_track('Stan'), 'The Marshall Mathers LP', 'could get cd by track';
    is_deeply [get_cds_by_artist('Michael Jackson')], ['Thriller', 'Bad'], 'could get cds by artist';
    
    is get_artist_by_track('Dirty Diana'), 'Michael Jackson', 'could get artist by track';
    is get_artist_by_cd('The Marshall Mathers LP'), 'Eminem', 'could get artist by cd';
    
    sub get_tracks_by_cd {
        my $cdtitle = shift;
        my @titles;
        my $rs = Schema->resultset('Track')->search(
            {
                'cd.title' => $cdtitle
            },
            {
                join     => [qw/ cd /],
            }
        );
        while (my $track = $rs->next) {
            push(@titles, $track->title);
        }
        return @titles;
    }
    
    sub get_tracks_by_artist {
        my $artistname = shift;
        my @titles;
        my $rs = Schema->resultset('Track')->search(
            {
                'artist.name' => $artistname
            },
            {
            join => {
                'cd' => 'artist'
            },
            }
        );
        while (my $track = $rs->next) {
            push(@titles, $track->title);
        }
        return @titles;
    }
    
    sub get_cd_by_track {
        my $tracktitle = shift;
        my $rs = Schema->resultset('CD')->search(
            {
                'tracks.title' => $tracktitle
            },
            {
                join     => [qw/ tracks /],
            }
        );
        my $cd = $rs->first;
        return $cd->title;
    }
    
    sub get_cds_by_artist {
        my $artistname = shift;
        my @titles;
        my $rs = Schema->resultset('CD')->search(
            {
                'artist.name' => $artistname
            },
            {
                join     => [qw/ artist /],
            }
        );
        while (my $cd = $rs->next) {
            push(@titles, $cd->title);
        }
        return @titles;
    }
    
    sub get_artist_by_track {
        my $tracktitle = shift;
        my $rs = Schema->resultset('Artist')->search(
            {
                'tracks.title' => $tracktitle
            },
            {
                join => {'cds' => 'tracks'}
            }
        );
        my $artist = $rs->first;
        return $artist->name;
    }
    
    sub get_artist_by_cd {
        my $cdtitle = shift;
        my $rs = Schema->resultset('Artist')->search(
            {
                'cds.title' => $cdtitle
            },
            {
                join     => [qw/ cds /],
            }
        );
        my $artist = $rs->first;
        return $artist->name;
    }
}

done_testing;
exit;
