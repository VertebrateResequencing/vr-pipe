#!/usr/bin/env perl

use strict;
use warnings;

my @db_settings = ("DBI:mysql:host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT};database=vertres_dbixclass;", $ENV{VRTRACK_RW_USER}, $ENV{VRTRACK_PASSWORD});

use VRPipe::Persistent::Schema;

my $schema = VRPipe::Persistent::Schema->connect(@db_settings);

my $deploy = 0;
$schema->deploy if $deploy;

my $insert = 0;
if ($insert) {
    my @artists = (['Michael Jackson'], ['Eminem']);
    $schema->populate('Artist', [
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
      my $artist = $schema->resultset('Artist')->find({
        name => $albums{$lp}
      });
      push @cds, [$lp, $artist->id];
    }
    
    $schema->populate('CD', [
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
      my $cdname = $schema->resultset('CD')->find({
        title => $tracks{$track},
      });
      push @tracks, [$cdname->id, $track];
    }
    
    $schema->populate('Track',[
      [qw/cd title/],
      @tracks,
    ]);
}

my $query = 0;
if ($query) {
    get_tracks_by_cd('Bad');
    get_tracks_by_artist('Michael Jackson');
    
    get_cd_by_track('Stan');
    get_cds_by_artist('Michael Jackson');
    
    get_artist_by_track('Dirty Diana');
    get_artist_by_cd('The Marshall Mathers LP');
    
    sub get_tracks_by_cd {
      my $cdtitle = shift;
      print "get_tracks_by_cd($cdtitle):\n";
      my $rs = $schema->resultset('Track')->search(
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
      my $rs = $schema->resultset('Track')->search(
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
      my $rs = $schema->resultset('CD')->search(
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
      my $rs = $schema->resultset('CD')->search(
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
      my $rs = $schema->resultset('Artist')->search(
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
      my $rs = $schema->resultset('Artist')->search(
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

my $constraint_checks = 1;
if ($constraint_checks) {
    print "calling create({name => 'Footon'})\n";
    my $artist = $schema->resultset('Artist')->create({name => 'Footon'});
    print "calling name('Barton')\n";
    $artist->name('Barton');
    my @cols = $artist->get_dirty_columns;
    print "dirty cols: (@cols)\n";
    $artist->update;
    @cols = $artist->get_dirty_columns;
    print "dirty cols: (@cols)\n";
    print "calling name()\n";
    print $artist->name, "\n";
    my @ids = $artist->id;
    print "ids: (@ids)\n";
    #$artist->delete;
    #use Data::Dumper;
    #print Dumper($artist), "\n";
}

exit;
