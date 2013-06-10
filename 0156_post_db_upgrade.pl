#!/usr/bin/env perl
use strict;
use warnings;
use JSON::XS;
use DBI;
use VRPipe;

warn "I will alter this database: ", VRPipe::Persistent::SchemaBase->get_dsn, " [ctrl-c now to abort]\n";
sleep(3);

my $json  = JSON::XS->new->utf8->canonical;
my $dbh   = DBI->connect(VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password, { RaiseError => 1, AutoCommit => 0 }) or die $DBI::errstr;
my $limit = 10000;

warn "\n ## First, I'll fill in File keyvallist column values based on the temporary file metadata table\n\n";
my $sth    = $dbh->prepare(q[update file set keyvallist = ? where id = ?]);
my $offset = 0;
my $total  = 0;
while (1) {
    my $select = $dbh->prepare(qq[SELECT f.file, f.metadata FROM ( SELECT file FROM temp_file_metadata LIMIT $offset, $limit ) o join temp_file_metadata f on f.file = o.file]);
    $select->execute;
    my $count = 0;
    # create all KeyValLists first, separately from updating the column, or we
    # get strange table lock issues
    my @execute_args;
    while (my $row = $select->fetch) {
        $total++;
        my ($id, $meta) = @$row;
        $meta = $json->decode($meta);
        my $kvl = VRPipe::KeyValList->get(hash => $meta);
        push(@execute_args, [$kvl->id, $id]);
        $count++;
    }
    warn " created/looked up $count keyvallists\n";
    
    $count = 0;
    eval {
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
    last if $count < $limit;
    
    $offset += $limit;
}

warn "\n ## Second, I'll check that worked and then drop the temporary file metadata table\n\n";
my $defined = VRPipe::File->search({ keyvallist => { '!=' => 0 } });
die "$defined file rows had a set keyvallist, yet $total files should have metadata?!\n" if $defined != $total;
my $drop_table = $dbh->prepare(q[DROP TABLE if exists temp_file_metadata]);
$drop_table->execute();
$dbh->commit;

warn "\n ## Third, I'll fill in DataElement keyvallist and filelist column values based on the temporary dataelement result table\n\n";
$sth    = $dbh->prepare(q[update dataelement set keyvallist = ?, filelist = ? where id = ?]);
$offset = 0;
while (1) {
    my $select = $dbh->prepare(qq[SELECT dataelement, metadata, files FROM temp_dataelement_result LIMIT $offset, $limit]);
    $select->execute;
    my $count = 0;
    my @execute_args;
    while (my $row = $select->fetch) {
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
    
    $count = 0;
    eval {
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
    last if $count < $limit;
    
    $offset += $limit;
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
