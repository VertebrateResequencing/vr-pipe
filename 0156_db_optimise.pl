#!/usr/bin/env perl
use strict;
use warnings;
use JSON::XS;
use Storable qw(thaw);
use DBI;
use VRPipe;

# we have for example over 40million File rows, but only 7million
# stepoutputfiles and 2million persistentarraymember files (probably datasource
# input files), and it might take days to do db_pre_upgrade if we don't first
# delete all the extraneous rows we've built up over time.

warn "I will alter this database: ", VRPipe::Persistent::SchemaBase->get_dsn, " [ctrl-c now to abort]\n";
sleep(3);

my $json  = JSON::XS->new->utf8->canonical;
my $dbh   = DBI->connect(VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password, { RaiseError => 1, AutoCommit => 0 }) or die $DBI::errstr;
my $limit = 10000;

warn "\n ## I'll delete File rows that aren't needed, even while the database is being used in production\n\n";
# first get the id of the most recent file row, so we won't delete any files
# created after we started this script
my $select_max = $dbh->prepare("SELECT id FROM file ORDER BY id DESC LIMIT 1");
$select_max->execute;
my $max_row = $select_max->fetch;
my ($max_id) = @$max_row;
warn "I won't delete any files with an id higher than $max_id (the current highest)\n";

my %needed_files;
warn "Will find datasource source files\n";
foreach my $ds (VRPipe::DataSource->search({})) {
    my $instance = $ds->_source_instance;
    next unless $instance->can('source_file');
    # skip bad source files
    my $source = $ds->source;
    next unless $source =~ /^\S+/;
    my $has_ps = VRPipe::PipelineSetup->search({ datasource => $ds->id });
    next unless $has_ps;
    $needed_files{ $instance->source_file->id } = 1;
}

warn "Will find datasource input files\n";
my $pager = VRPipe::PersistentArrayMember->get_column_values_paged('class_id', { class => 'VRPipe::File' });
while (my $fids = $pager->next) {
    foreach my $fid (@$fids) {
        $needed_files{$fid} = 1;
    }
}

warn "Will find step output files\n";
$pager = VRPipe::StepOutputFile->get_column_values_paged('file', {}); # we can't skip temp files or this table will reference non-existent files
while (my $fids = $pager->next) {
    foreach my $fid (@$fids) {
        $needed_files{$fid} = 1;
    }
}

my $num_needed = keys %needed_files;
warn "Will find the chain of files leading to or from our needed files ($num_needed)\n";
my %chained_files;
my $chained = 0;
while (my ($fid) = each %needed_files) {
    add_chained_ids($fid);
    $chained++;
    if ($chained % $limit == 0) {
        warn " found chain for $chained files so far\n";
    }
}

sub add_chained_ids {
    my $fid = shift;
    my $pm = VRPipe::File->get_column_values(['parent', 'moved_to'], { id => $fid }, { rows => 1 });
    my ($parent, $moved_to) = @{ $pm->[0] };
    if ($parent) {
        $chained_files{$parent} = 1;
        add_chained_ids($parent);
    }
    if ($moved_to) {
        $chained_files{$moved_to} = 1;
        while (1) {
            ($moved_to) = VRPipe::File->get_column_values('moved_to', { id => $moved_to }, { rows => 1 });
            $moved_to || last;
            $chained_files{$moved_to} = 1;
        }
    }
    
    my $current_id = $fid;
    while (1) {
        my ($moved_from) = VRPipe::File->get_column_values('id', { moved_to => $current_id }, { rows => 1 });
        $moved_from || last;
        $chained_files{$moved_from} = 1;
        $current_id = $moved_from;
    }
}

$dbh = DBI->connect(VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password, { RaiseError => 1, AutoCommit => 0 }) or die $DBI::errstr;
$dbh->do(q[SET FOREIGN_KEY_CHECKS = 0]);
my $sth = $dbh->prepare("delete from file where id = ?");
warn "Will delete File rows no longer needed\n";
my $offset = 0;
my $last_id;
while (1) {
    # (order by id could be used to help avoid hitting foreign key contraints,
    # but it never worked correctly anyway, and it results in >500second
    # selects)
    my $select = $dbh->prepare("SELECT id, metadata FROM file where id < $max_id LIMIT $offset, $limit");
    $select->execute;
    
    my $deleted = 0;
    my $skipped = 0;
    my $count   = 0;
    eval {
        while (my $row = $select->fetch) {
            $count++;
            
            my ($id, $meta) = @$row;
            if (exists $needed_files{$id} || exists $chained_files{$id}) {
                $skipped++;
                next;
            }
            
            $meta = thaw($meta);
            if (keys %{$meta}) {
                $skipped++;
                next;
            }
            
            $last_id = $id;
            $sth->execute($id);
            $deleted++;
        }
        $dbh->commit;
        warn " deleted $deleted rows (kept $skipped; offset $offset)\n";
    };
    if ($@) {
        my $error = $@;
        eval { $dbh->rollback };
        die "Transaction aborted at File row with id $last_id because $error";
    }
    last if $count < $limit;
    
    $offset += $limit;
}

$dbh->do(q[SET FOREIGN_KEY_CHECKS = 1]);
$dbh->disconnect;
warn "\n ## All done; now turn off vrpipe-server, kill any remaining jobs and run 0156_pre_db_upgrade\n\n";
exit;
