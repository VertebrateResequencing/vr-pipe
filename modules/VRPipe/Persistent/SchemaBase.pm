
=head1 NAME

VRPipe::Persistent::SchemaBase - the backend for connecting to the correct     
                            database used for Persistent objects

=head1 SYNOPSIS
    
    use VRPipe::Base;
    
    class VRPipe::Persistent::Schema extends VRPipe::Persistent::SchemaBase {
        __PACKAGE__->load_namespaces(result_namespace =>
                                        ['+VRPipe::DirA', '+VRPipe::DirB']);
    }

=head1 DESCRIPTION

This is a subclass of DBIx::Class::Schema that can default connection details
to site-wide configuration values (from L<VRPipe::Config>).

To choose between using the production or testing database that has been
configured, call
C<VRPipe::Persistent::SchemaBase->database_deployment('testing')> prior to
calling C<connect()> on a subclass of this class.

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

class VRPipe::Persistent::SchemaBase extends (DBIx::Class::Schema, VRPipe::Base::Moose) {
    use MooseX::NonMoose;
    use VRPipe::Config;
    use File::HomeDir;
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
            $dsn  = $thing->{dsn};
            $user = $thing->{user} || $class->get_user;
            $pass = $thing->{password} || $class->get_password;
        }
        else {
            $dsn  = $thing || $class->get_dsn;
            $user = shift  || $class->get_user;
            $pass = shift  || $class->get_password;
        }
        
        my $args = shift || { AutoCommit => 1, RaiseError => 1, PrintError => 0 };
        
        my $dbtype = $class->get_dbtype;
        $args->{sqlite_use_immediate_transaction} = 1 if $dbtype =~ /sqlite/i;
        
        $class->SUPER::connect($dsn, $user, $pass, $args);
    }
    
    method database_deployment (ClassName $class: Str $set?) {
        if ($set) {
            if ($set eq 'production' || $set eq 'testing') {
                $DATABASE_DEPLOYMENT = $set;
            }
            else {
                die "Invalid deployment '$set'\n";
            }
        }
        return $DATABASE_DEPLOYMENT;
    }
    
    method get_dsn (ClassName $class:) {
        my %details;
        foreach my $detail (qw(dbtype dbname dbhost dbport)) {
            my $method_name = $DATABASE_DEPLOYMENT . '_' . $detail;
            $details{$detail} = $vrp_config->$method_name();
        }
        
        my $dbtype = $details{dbtype};
        if ($dbtype =~ /mysql/i) {
            unless (defined $details{dbname} && defined $details{dbhost} && defined $details{dbport}) {
                die "SiteConfig did not have all the necessary database details specified; did modules/VRPipe/SiteConfig.pm get created and installed properly?\n";
            }
            return "dbi:mysql:database=$details{dbname};host=$details{dbhost};port=$details{dbport}";
        }
        elsif ($dbtype =~ /pg|postgres/i) {
            unless (defined $details{dbname} && defined $details{dbhost} && defined $details{dbport}) {
                die "SiteConfig did not have all the necessary database details specified; did modules/VRPipe/SiteConfig.pm get created and installed properly?\n";
            }
            return "dbi:Pg:dbname=$details{dbname};host=$details{dbhost};port=$details{dbport}";
        }
        elsif ($dbtype =~ /sqlite/i) {
            unless (defined $details{dbname}) {
                die "SiteConfig did not have all the necessary database details specified; did modules/VRPipe/SiteConfig.pm get created and installed properly?\n";
            }
            my $name = $details{dbname};
            if ($name =~ /^~/) {
                my $home;
                if ($name =~ /^~([\w]+)/) {
                    $home = File::HomeDir->users_home($1);
                }
                else {
                    $home = File::HomeDir->my_home;
                }
                $name =~ s/^~[\w]*/$home/;
            }
            return "dbi:SQLite:dbname=$name";
        }
        else {
            die "Unsupported database type '$dbtype'\n";
        }
    }
    
    method get_dbtype (ClassName $class:) {
        my $method_name = $DATABASE_DEPLOYMENT . '_dbtype';
        my $dbtype      = $vrp_config->$method_name();
        $dbtype = "$dbtype";
        
        # correct capitlisation for the benefit of DBIx::Class::DeploymentHandler
        if ($dbtype =~ /mysql/i) {
            $dbtype = 'MySQL';
        }
        elsif ($dbtype =~ /sqlite/i) {
            $dbtype = 'SQLite';
        }
        
        return $dbtype;
    }
    
    method get_schema_dir (ClassName $class:) {
        my $dir = $vrp_config->schema_directory();
        return "$dir";
    }
    
    method get_user (ClassName $class:) {
        my $method_name = $DATABASE_DEPLOYMENT . '_username';
        my $user        = $vrp_config->$method_name();
        return "$user";  # incase $user is an Env object, force the stringification
    }
    
    method get_password (ClassName $class:) {
        my $method_name = $DATABASE_DEPLOYMENT . '_password';
        my $pass        = $vrp_config->$method_name();
        return "$pass";  # incase $pass is an Env object, force the stringification
    }
}

1;
