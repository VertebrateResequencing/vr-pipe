=head1 NAME

VRPipe::Persistent::SchemaBase - the backend for connecting to the correct
                                 database used for Persistent objects

=head1 SYNOPSIS

use VRPipe::Base;

class VRPipe::Persistent::Schema extends VRPipe::Persistent::SchemaBase {
    __PACKAGE__->load_namespaces(result_namespace => ['+VRPipe::DirA', '+VRPipe::DirB']);
}

=head1 DESCRIPTION

This is a subclass of DBIx::Class::Schema that can default connection details to
site-wide configuration values (from VRPipe::Config).

To choose between using the production or testing database that has been
configured, call VRPipe::Persistent::SchemaBase->database_deployment('testing')
prior to calling connect() on a subclass of this class.

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use VRPipe::Base;

class VRPipe::Persistent::SchemaBase extends (DBIx::Class::Schema, VRPipe::Base::Moose) {
    use MooseX::NonMoose;
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    
    our $DATABASE_DEPLOYMENT = 'production';
    
    # __PACKAGE__->exception_action(sub { My::ExceptionClass->throw(@_) });
    __PACKAGE__->stacktrace(1);
    
    sub new {
        return shift->connect(@_);
    }
    
    sub connect {
        my $class = shift;
        
        my ($dsn, $user, $pass);
        
        my $thing = shift;
        if (ref $thing) {
            $dsn = $thing->{dsn};
            $user = $thing->{user} || $class->get_user;
            $pass = $thing->{password} || $class->get_password;
        }
        else {
            $dsn = $thing || $class->get_dsn;
            $user = shift || $class->get_user;
            $pass = shift || $class->get_password;
        }
        
        my $args = shift || { AutoCommit => 1, RaiseError => 1, PrintError => 0 };
        
        $class->SUPER::connect($dsn, $user, $pass, $args);
    }
    
    method database_deployment (ClassName $class: Str $set?) {
        if ($set && ($set eq 'production' || $set eq 'testing')) {
            $DATABASE_DEPLOYMENT = $set;
        }
        return $DATABASE_DEPLOYMENT;
    }
    
    method get_dsn (ClassName $class:) {
        my %details;
        foreach my $detail (qw(dbtype dbname dbhost dbport)) {
            my $method_name = $DATABASE_DEPLOYMENT.'_'.$detail;
            $details{$detail} = $vrp_config->$method_name();
        }
        
        my $dbtype = $details{dbtype};
        if ($dbtype =~ /mysql/i) {
            return "dbi:mysql:database=$details{dbname};host=$details{dbhost};port=$details{dbport}";
        }
        elsif ($dbtype =~ /pg|postgres/i) {
            return "dbi:Pg:dbname=$details{dbname};host=$details{dbhost};port=$details{dbport}";
        }
        elsif ($dbtype =~ /sqlite/i) {
            return "dbi:SQLite:dbname=$details{dbname}";
        }
        else {
            die "Unsupported database type '$dbtype'\n";
        }
    }
    
    method get_user (ClassName $class:) {
        my $method_name = $DATABASE_DEPLOYMENT.'_username';
        my $user = $vrp_config->$method_name();
        return "$user"; # incase $user is an Env object, force the stringification
    }
    
    method get_password (ClassName $class:) {
        my $method_name = $DATABASE_DEPLOYMENT.'_password';
        my $pass = $vrp_config->$method_name();
        return "$pass"; # incase $pass is an Env object, force the stringification
    }
}

1;
