#!/usr/bin/env perl
use strict;
use warnings;
use JSON::XS;
use DBI;
use VRPipe;

warn "I will alter this database: ", VRPipe::Persistent::SchemaBase->get_dsn, " [ctrl-c now to abort]\n";
sleep(3);

my $json = JSON::XS->new->utf8->canonical;
my $dbh = DBI->connect(VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password, { RaiseError => 1, AutoCommit => 0 }) or die $DBI::errstr;

warn "\n ## First, I'll fill in File keyvallist column values based on the temporary file metadata table\n\n";
my $select = $dbh->prepare(q[SELECT file, metadata FROM temp_file_metadata]);
$select->execute;
my $sth   = $dbh->prepare(q[update file set keyvallist = ? where id = ?]);
my $total = 0;
while (1) {
    # create all KeyValLists first, separately from updating the column, or we
    # get strange table lock issues
    my $last  = 0;
    my $count = 0;
    my @execute_args;
    for (1 .. 10) {
        my $row = $select->fetch;
        unless ($row) { $last = 1; last; }
        $total++;
        my ($id, $meta) = @$row;
        $meta = $json->decode($meta);
        my $kvl = VRPipe::KeyValList->get(hash => $meta);
        push(@execute_args, [$kvl->id, $id]);
        $count++;
    }
    warn " created/looked up $count keyvallists\n";
    
    eval {
        $count = 0;
        foreach (@execute_args) {
            $sth->execute(@$_);
            $count++;
        }
        $dbh->commit;
        warn " filled in $count rows\n";
    };
    if ($@) {
        my $error = $@;
        eval { $dbh->rollback };
        die "Transaction aborted because $error";
    }
    
    last if $last;
}

warn "\n ## Second, I'll check that worked and then drop the temporary file metadata table\n\n";
my $defined = VRPipe::File->search({ keyvallist => { '!=' => 0 } });
die "$defined file rows had a set keyvallist, yet $total files should have metadata?!\n" if $defined != $total;
my $drop_table = $dbh->prepare(q[DROP TABLE if exists temp_file_metadata]);
$drop_table->execute();
$dbh->commit;

warn "\n ## Third, I'll fill in DataElement keyvallist and filelist column values based on the temporary dataelement result table\n\n";
$select = $dbh->prepare(q[SELECT dataelement, metadata, files FROM temp_dataelement_result]);
$select->execute;
$sth = $dbh->prepare(q[update dataelement set keyvallist = ?, filelist = ? where id = ?]);
while (1) {
    my $last  = 0;
    my $count = 0;
    my @execute_args;
    for (1 .. 10) {
        my $row = $select->fetch;
        unless ($row) { $last = 1; last; }
        my ($id, $meta, $files) = @$row;
        $meta = $json->decode($meta);
        my $kvl = VRPipe::KeyValList->get(hash => $meta);
        $files = $json->decode($files);
        my @files;
        foreach my $fid (@$files) {
            push(@files, VRPipe::File->get(id => $fid));
        }
        my $fl = VRPipe::FileList->create(files => \@files);
        push(@execute_args, [$kvl->id, $fl->id, $id]);
        $count++;
    }
    warn " created/looked up $count keyvallists and filelists\n";
    
    eval {
        $count = 0;
        foreach (@execute_args) {
            $sth->execute(@$_);
            $count++;
        }
        $dbh->commit;
        warn " filled in $count rows\n";
    };
    if ($@) {
        my $error = $@;
        eval { $dbh->rollback };
        die "Transaction aborted because $error";
    }
    
    last if $last;
}
warn "\n ## Finally, I'll check that worked and then drop the temporary dataelement result table\n\n";
my $undefs = VRPipe::DataElement->search({ keyvallist => 0 });
die "$undefs dataelement rows still had a keyvallist of 0?!\n" if $undefs;
$undefs = VRPipe::DataElement->search({ filelist => 0 });
die "$undefs dataelement rows still had a filelist of 0?!\n" if $undefs;
$drop_table = $dbh->prepare(q[DROP TABLE if exists temp_dataelement_result]);
$drop_table->execute();
$dbh->commit;
$dbh->disconnect;

warn "## All done; you can now use VRPipe as normal\n\n";
exit;
