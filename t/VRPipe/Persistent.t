#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 9;
    
    use_ok('VRPipe::Persistent::Schema');
    use_ok('VRPipe::Persistent');
    use_ok('File::Spec');
    
    #my @db_settings = ("DBI:mysql:host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT};database=vertres_dbixclass;", $ENV{VRTRACK_RW_USER}, $ENV{VRTRACK_PASSWORD});
    use Test::DBIx::Class {
        schema_class => 'VRPipe::Persistent::Schema',
        connect_info => ['dbi:SQLite:dbname=:memory:','','']
    }; #, 'VRPipe::Artist' => {-as => 'Artist'}, 'VRPipe::CD' => {-as => 'CD'}, 'VRPipe::Track' => {-as => 'Track'}
}

# some quick basic tests
ok my $resultset = ResultSet('Artist'), "Yeah, some Artists in the database";

ok my $artist = $resultset->create({name => 'Barton'}), 'Created Barton';
ok $artist = $resultset->create({name => 'Footon'}), 'Created Footon';

ok my $footon = $resultset->search({'name' => 'Footon'})->first, 'Searched for Footon';

is_fields [qw/artistid name/], $footon, [2, 'Footon'], 'Footon has the expected fields';


# was playing with some contraint checking
my $constraint_checks = 0;
if ($constraint_checks) {
    print "calling name('Barton')\n";
    use TryCatch;
    try {
        $artist->name('Barton');
    }
    catch (MooseX::Error::Exception::Class $error) {
        warn "got error $error\n";
    }
    catch {
        die "unhandled error\n";
    }
    
    my @cols = $artist->get_dirty_columns;
    print "dirty cols: (@cols)\n";
    $artist->update;
    @cols = $artist->get_dirty_columns;
    print "dirty cols: (@cols)\n";
    print "calling name()\n";
    print $artist->name, "\n";
    my @ids = $artist->id;
    print "ids: (@ids)\n";
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
