
=head1 NAME

VRPipe::Persistent::Schema - the frontend for connecting to and getting at     
                        persistently stored objects

=head1 SYNOPSIS
    
    use VRPipe::Persistent::Schema;
    
    my $schema = VRPipe::Persistent::Schema->connect();
    
    my $resultset = $schema->resultset('Job');
    
    # $resultset isa DBIx::Class::ResultSet; do stuff with it eg:
    
    while (my $job = $resultset->next) {
	print $job->id;
    }
    
    $finished_jobs_rs = $schema->resultset('Job')->search({ finished => 1 });

=head1 DESCRIPTION

This is a subclass of DBIx::Class::Schema that loads VRPipe Persistent classes
and lets you search for/create particular instances. The details of which
database to connect to are automatically determined from the site-wide
configuration.

See
L<http://search.cpan.org/~abraxxa/DBIx-Class-0.08127/lib/DBIx/Class/ResultSet.pm>.

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::Persistent::Schema extends VRPipe::Persistent::SchemaBase {
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::Persistent::ConverterFactory;
    
    our $VERSION = 26;
    __PACKAGE__->load_classes({
            'VRPipe' => [
                qw/Step Scheduler Job Requirements
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
                  StepStats FarmServer Runner MessageTracker/
            ]
        }
    );
    
    # deploy method overridden in order to add indexes in a db-dependent manner
    sub deploy {
        my ($self, $sqltargs, $dir) = @_;
        
        $self->throw_exception("Can't deploy without storage") unless $self->storage;
        $self->storage->deploy($self, undef, $sqltargs, $dir);
        
        my $idx_cmds = $self->get_idx_sql('create');
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
        my ($self, $mode) = @_;
        $mode = 'create' unless $mode;
        
        my $dbtype = lc(VRPipe::Persistent::SchemaBase->get_dbtype);
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {});
        
        return $converter->index_statements($self, $mode);
    }
    
    sub get_db_schema_version {
        my ($self) = @_;
        my $dbtype = lc(VRPipe::Persistent::SchemaBase->get_dbtype);
        my $converter = VRPipe::Persistent::ConverterFactory->create($dbtype, {});
        return $converter->get_db_schema_version($self);
    }
}

1;
