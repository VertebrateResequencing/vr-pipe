#!/usr/bin/env perl

use strict;
use warnings;
use Tangram::Relational;
use DBI;
use VRPipe::Tangram;
use Data::Dumper;

my @db_settings = ("DBI:mysql:host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT};database=vertres_tangram;", $ENV{VRTRACK_RW_USER}, $ENV{VRTRACK_PASSWORD});

if (0) {
    print "connecting...\n";
    my $dbh = DBI->connect(@db_settings);
    print "deploying...\n";
    Tangram::Relational->deploy($VRPipe::Tangram::schema, $dbh);
}

print "acccessing storage...\n";
my $storage = VRPipe::Tangram::db(@db_settings);
print Dumper($storage);
#print "attempting to connect to storage...\n";
#$storage->connect;

print "creating object...\n";
my $object = bless { first_name => "Homer", last_lame => "Simpson" }, "NaturalPerson";
print Dumper($object);
print "inserting object...\n";
my $oid = $storage->insert($object);
# *** failing at insert, or at connect above

print "loading object...\n";
my $homer = $storage->load($oid);
print Dumper($homer);


exit;
