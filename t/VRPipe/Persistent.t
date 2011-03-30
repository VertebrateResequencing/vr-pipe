#!/usr/bin/env perl
use strict;
use warnings;
use TryCatch;


BEGIN {
    use Test::Most tests => 18;
    
    use_ok('VRPipe::Persistent::Schema');
    use_ok('VRPipe::Persistent');
    use_ok('File::Spec');
    
    use_ok('VRPipe::Test');
    
    use Test::DBIx::Class {
        schema_class => 'VRPipe::Persistent::Schema',
        force_drop_table => 1,
        keep_db => 1,
        connect_info => ["dbi:mysql:host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT};database=vertres_dbixclass;", $ENV{VRTRACK_RW_USER}, $ENV{VRTRACK_PASSWORD}] #['dbi:SQLite:dbname=:memory:','','']
    }; #, 'VRPipe::Artist' => {-as => 'Artist'}, 'VRPipe::CD' => {-as => 'CD'}, 'VRPipe::Track' => {-as => 'Track'}
}

my $test = VRPipe::Test->new;
$test->first_name('foo');
is $test->first_name, 'foo', 'Made a VRPipe::Test and it worked';
is $test->foo, 'food', 'base foo method worked';
is $test->verbose, 0, 'verbose defaults to 0';
$test->warn("debuggable warning");
$test->verbose(1);
$test->warn("debuggable warning at verbose 1");

# some quick basic tests
ok my $resultset = ResultSet('Artist'), "Yeah, some Artists in the database";

ok my $artist = $resultset->create({name => 'Barton'}), 'Created Barton';
ok $artist = $resultset->create({name => 'Footon'}), 'Created Footon';

$artist->warn("debuggable warning");
$artist->verbose(1);
$artist->warn("debuggable warning at verbose 1");
is $artist->foo, 'food', 'base foo method worked on an artist';

ok my $footon = $resultset->search({'name' => 'Footon'})->first, 'Searched for Footon';

is_fields [qw/artistid name/], $footon, [2, 'Footon'], 'Footon has the expected fields';


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
try {
    $artist->name('12345678901234567890123456789012345678901234567890123456789012345');
}
catch (MooseX::Error::Exception::Class $error) {
    warn "got error $error\n";
}
catch ($error) {
    die "unhandled error $error\n";
}


# relationship testing from dbic example code
reset_schema; # this counts as a test
{
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
}

{
    get_tracks_by_cd('Bad');
    get_tracks_by_artist('Michael Jackson');
    
    get_cd_by_track('Stan');
    get_cds_by_artist('Michael Jackson');
    
    get_artist_by_track('Dirty Diana');
    get_artist_by_cd('The Marshall Mathers LP');
    
    sub get_tracks_by_cd {
      my $cdtitle = shift;
      print "get_tracks_by_cd($cdtitle):\n";
      my $rs = Schema->resultset('Track')->search(
        {
          'cd.title' => $cdtitle
        },
        {
          join     => [qw/ cd /],
        }
      );
      while (my $track = $rs->next) {
        print $track->title . "\n";
      }
      print "\n";
    }
    
    sub get_tracks_by_artist {
      my $artistname = shift;
      print "get_tracks_by_artist($artistname):\n";
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
        print $track->title . "\n";
      }
      print "\n";
    }
    
    sub get_cd_by_track {
      my $tracktitle = shift;
      print "get_cd_by_track($tracktitle):\n";
      my $rs = Schema->resultset('CD')->search(
        {
          'tracks.title' => $tracktitle
        },
        {
          join     => [qw/ tracks /],
        }
      );
      my $cd = $rs->first;
      print $cd->title . "\n\n";
    }
    
    sub get_cds_by_artist {
      my $artistname = shift;
      print "get_cds_by_artist($artistname):\n";
      my $rs = Schema->resultset('CD')->search(
        {
          'artist.name' => $artistname
        },
        {
          join     => [qw/ artist /],
        }
      );
      while (my $cd = $rs->next) {
        print $cd->title . "\n";
      }
      print "\n";
    }
    
    sub get_artist_by_track {
      my $tracktitle = shift;
      print "get_artist_by_track($tracktitle):\n";
      my $rs = Schema->resultset('Artist')->search(
        {
          'tracks.title' => $tracktitle
        },
        {
          join => {
            'cds' => 'tracks'
          }
        }
      );
      my $artist = $rs->first;
      print $artist->name . "\n\n";
    }
    
    sub get_artist_by_cd {
      my $cdtitle = shift;
      print "get_artist_by_cd($cdtitle):\n";
      my $rs = Schema->resultset('Artist')->search(
        {
          'cds.title' => $cdtitle
        },
        {
          join     => [qw/ cds /],
        }
      );
      my $artist = $rs->first;
      print $artist->name . "\n\n";
    }
}


done_testing;
exit;
