=head1 NAME

VRPipe::Persistent::Schema - the frontend for connecting to and getting at
                             persistently stored objects

=head1 SYNOPSIS

use VRPipe::Persistent::Schema;

my $schema = VRPipe::Persistent::Schema->connect();

my $resultset = $schema->resultset('Job');

# $resultset isa DBIx::Class::ResultSet; do stuff with it
# http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSet.pm
# eg:

while (my $job = $resultset->next) {
    print $job->id;
}

$finished_jobs_resultset = $schema->resultset('Job')->search({ finished => 1 });

=head1 DESCRIPTION

This is a subclass of DBIx::Class::Schema that loads VRPipe Persistent classes
and lets you search for/create particular instances. The details of which
database to connect to are automatically determined from the site-wide
configuration.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

class VRPipe::Persistent::Schema extends VRPipe::Persistent::SchemaBase {

    use VRPipe::Persistent::SchemaBase;
    use VRPipe::Persistent::ConverterFactory;

    our $VERSION = 19;
    __PACKAGE__->load_classes({'VRPipe' => [qw/Step Scheduler Job Requirements
                                               DataSource DataElement Pipeline
                                               StepCmdSummary StepMember File
                                               PipelineSetup StepBehaviour
                                               StepState Submission StepAdaptor
                                               PersistentArray StepOption
                                               PersistentArrayMember Manager
                                               StepIODefinition StepOutputFile
                                               DataElementState DataElementLink
                                               LocalSchedulerJob
                                               LocalSchedulerJobState
                                               StepStats/]});

    # deploy method overridden in order to add indexes in a db-dependent manner
    sub deploy {
        my ($self, $sqltargs, $dir) = @_;

        $self->throw_exception("Can't deploy without storage") unless $self->storage;
        $self->storage->deploy($self, undef, $sqltargs, $dir);

        my $dbtype = lc(VRPipe::Persistent::SchemaBase->get_dbtype);
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {});

        my $idx_cmds = $converter->index_statements($self,'create');

        if ($idx_cmds) {
            $self->storage->dbh_do(
                sub {
                    my ($storage, $dbh, $idx_cmds) = @_;
                    foreach my $cmd (@{$idx_cmds}) {
                       $dbh->do($cmd);
                    }
                },
                $idx_cmds
            );
        }
    }

    sub get_idx_sql {
        my ($self,$mode) = @_;
		$mode = 'create' unless $mode;

        my $dbtype = lc(VRPipe::Persistent::SchemaBase->get_dbtype);
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {});

##      my $idx_cmds = $converter->index_statements($self, $mode);	# Does not work, self is just a string 'VRPipe::Persistent::Schema'

		my $idx_cmds;
        foreach my $class (keys %{$self->class_mappings}) { # self is a hashref with class_mappings and meta data !!
            my $table_name = $class;
            $table_name =~ s/.*:://;
            $table_name = lc($table_name);
            my $meta = $class->meta;
            my $for_indexing = $meta->get_attribute('idx_keys')->get_value($meta);

            if (keys %{$for_indexing}) {
				push(@{$idx_cmds}, @{$converter->get_index_statements($table_name, $for_indexing, $mode)});
            }
        }
		return $idx_cmds;
    }

}

1;
