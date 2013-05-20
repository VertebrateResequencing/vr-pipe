#!/usr/bin/env perl
use strict;
use warnings;
use JSON::XS;
use Storable qw(thaw);
use DBI;
use VRPipe;

warn "I will alter this database: ", VRPipe::Persistent::SchemaBase->get_dsn, " [ctrl-c now to abort]\n";
sleep(3);

my $json = JSON::XS->new->utf8->canonical;
my $dbh = DBI->connect(VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password, { RaiseError => 1, AutoCommit => 0 }) or die $DBI::errstr;

warn "\n ## First, I'll do in-place updates of columns containing Storable values, converting them to JSON values\n\n";
foreach ([qw(Job output_files)], [qw(Step options_definition inputs_definition outputs_definition)], [qw(DataSource options)], [qw(PipelineSetup options)], [qw(Requirements custom)], [qw(StepAdaptor adaptor_hash)], [qw(StepIODefinition metadata)]) {
    my ($table, @cols) = @$_;
    my $class = "VRPipe::$table";
    $table = lc($table);
    my $select = $dbh->prepare("SELECT id, " . join(", ", @cols) . " FROM $table");
    $select->execute;
    my $sth = $dbh->prepare("update $table set " . join(", ", map { "`$_` = ?" } @cols) . " where id = ?");
    warn "Will update table $table, columns ", join(", ", @cols), "\n";
    while (1) {
        my $last = 0;
        eval {
            my $count = 0;
            for (1 .. 1000) {
                my $row = $select->fetch;
                unless ($row) { $last = 1; last; }
                my ($id, @orig_vals) = @$row;
                my @new_vals;
                foreach my $val (@orig_vals) {
                    eval { $val = thaw($val); };
                    if ($@) {
                        my $error = $@;
                        if ($error =~ /Storable binary image .+ more recent than/) {
                            # we already converted this on a previous try?
                            next;
                        }
                        die $error;
                    }
                    push(@new_vals, $json->encode($val));
                }
                next if @new_vals == 0;
                if (@new_vals != @cols) { die "Looks like $table row with id $id is in a bad state; fix it manually before trying again\n"; }
                $sth->execute(@new_vals, $id);
                $count++;
            }
            $dbh->commit;
            warn " updated $count rows\n";
        };
        if ($@) {
            my $error = $@;
            eval { $dbh->rollback };
            die "Transaction aborted because $error";
        }
        last if $last;
    }
}

warn "\n ## Second, I'll place current File metadata column values into a temporary table\n\n";
my $create_table = $dbh->prepare(q[CREATE TABLE if not exists temp_file_metadata (file INT NOT NULL, PRIMARY KEY (file), metadata TEXT)]);
$create_table->execute();
warn "Created table temp_file_metadata; will store metadata from file table in it...\n";
my $select = $dbh->prepare(q[SELECT id, metadata FROM file]);
$select->execute;
my $sth = $dbh->prepare(q[INSERT INTO temp_file_metadata (file, metadata) VALUES (?, ?) ON DUPLICATE KEY UPDATE metadata = ?]);
while (1) {
    my $last = 0;
    eval {
        my $count = 0;
        for (1 .. 1000) {
            my $row = $select->fetch;
            unless ($row) { $last = 1; last; }
            my ($id, $meta) = @$row;
            $meta = thaw($meta);
            delete $meta->{original_pg_chain}; # it is large and unnecessary to keep
            $meta = $json->encode($meta);
            $sth->execute($id, $meta, $meta);
            $count++;
        }
        $dbh->commit;
        warn " stored $count rows\n";
    };
    if ($@) {
        my $error = $@;
        eval { $dbh->rollback };
        die "Transaction aborted because $error";
    }
    last if $last;
}

warn "\n ## Finally, I'll place current DataElement result column values into a temporary table\n\n";
$create_table = $dbh->prepare(q[CREATE TABLE if not exists temp_dataelement_result (dataelement INT NOT NULL, PRIMARY KEY (dataelement), metadata TEXT, files TEXT)]);
$create_table->execute();
warn "Created table temp_dataelement_result; will store result metadata and input files from dataelement table in it...\n";
$select = $dbh->prepare(q[SELECT id, result FROM dataelement]);
$select->execute;
$sth = $dbh->prepare(q[INSERT INTO temp_dataelement_result (dataelement, metadata, files) VALUES (?, ?, ?) ON DUPLICATE KEY UPDATE metadata = ?, files = ?]);
my $get_files_h = $dbh->prepare(q[SELECT class_id from persistentarraymember where persistentarray = ?]);
while (1) {
    my $last = 0;
    eval {
        my $count = 0;
        for (1 .. 1000) {
            my $row = $select->fetch;
            unless ($row) { $last = 1; last; }
            my ($id, $result) = @$row;
            $result = thaw($result);
            my $peristentarray = delete $result->{paths} || {};
            my @files;
            $get_files_h->execute($peristentarray);
            while (my $file_row = $get_files_h->fetch) {
                push(@files, $file_row->[0]);
            }
            my $meta  = $json->encode($result);
            my $files = $json->encode(\@files);
            $sth->execute($id, $meta, $files, $meta, $files);
            $count++;
        }
        $dbh->commit;
        warn " stored $count rows\n";
    };
    if ($@) {
        my $error = $@;
        eval { $dbh->rollback };
        die "Transaction aborted because $error";
    }
    last if $last;
}

warn "\n ## All done; now run vrpipe-db_upgrade followed by 0156_post_db_upgrade\n\n";
exit;
