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
    our $VERSION = 6;
    __PACKAGE__->load_classes({'VRPipe' => [qw/Step Scheduler Job Requirements
                                               DataSource DataElement Pipeline
                                               StepCmdSummary StepMember
                                               PipelineSetup StepBehaviour
                                               StepState Submission StepAdaptor
                                               PersistentArray StepOption
                                               PersistentArrayMember Manager
                                               File StepIODefinition/]});
}

1;
