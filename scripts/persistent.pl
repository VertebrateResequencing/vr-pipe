#!/usr/bin/env perl

use strict;
use warnings;

my @db_settings = ("DBI:mysql:host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT};database=vertres_persistent;", $ENV{VRTRACK_RW_USER}, $ENV{VRTRACK_PASSWORD});

use VRPipe::Persistent;

eval {
    #my $vrp = new VRPipe::Persistent('vrp.dbm', undef, 'DB_File');
    my $vrp = new VRPipe::Persistent(@db_settings, 'test');
    
    if ($vrp->restore('Fred', 'Flintstone')) {
        warn "Fred Flintstone was restored\n";
        if ($vrp->restore('Fred', 'Rubble')) {
            warn "Fred Rubble already exists:\n";
            $vrp->print;
        }
        else {
            warn "Will rename Flintstone to Rubble\n";
            $vrp->lastname('Rubble');
            $vrp->update();
            $vrp->print;
        }
    }
    else {
        warn "Fred Flintstone was not found, will create him\n";
        $vrp->firstname('Fred');
        $vrp->lastname('Flintstone');
        $vrp->telnum('650-555-1111');
        $vrp->age(45);
        $vrp->bday('1954-01-23 22:09:54');
        $vrp->save; # save behaving like insert; you can't do this twice in a row
        $vrp->print;
    }
};

if ($@) {
    warn "Exception caught: $@\n";
}


exit;
