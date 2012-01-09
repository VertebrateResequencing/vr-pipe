package VRPipeTest;
use strict;
use warnings;
use base 'Test::DBIx::Class';
use VRPipe::Persistent::SchemaBase;
use File::Spec;
use File::Which qw(which);

BEGIN { unshift(@INC, './modules') }
VRPipe::Persistent::SchemaBase->database_deployment('testing'); # do not change this under any circumstance!

use constant WIN32 => $^O eq 'MSWin32';
my $quote = WIN32 ? q/"/ : q/'/;

sub import {
    my $class  = shift;
    my $caller = caller(0);
    my %args = @_;
    
    # skip if required_env or required_exe not satisfied
    while (my ($key, $val) = each %args) {
        if ($key =~ /^required_(env|exe)/) {
            my $type = $1;
            my @to_check = ref($val) ? @{$val} : ($val);
            foreach (@to_check) {
                skip_all_unless_exists($type => $_);
            }
        }
    }
    
    # copy/pasted  Test::DBIx::Class to get the exports
    my ($schema_manager, $merged_config, @exports) = $class->_initialize();
    my $exporter = Sub::Exporter::build_exporter({
        exports => [
            Schema => sub {
                return sub {
                    return $schema_manager->schema;
                }
            }
        ],
        into_level => 1,    
    });
    
    $class->$exporter(qw/Schema/, @exports);
}

sub _initialize_schema {
    my $class = shift;
    
    my $config  = { schema_class => 'VRPipe::Persistent::Schema',
                    force_drop_table => 1,
                    keep_db => 1,
                    connect_info => [VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password] };
    
    return $class->SUPER::_initialize_schema($config);
}


# following methods were stolen and slightly modified from
# Test::Skip::UnlessExistsExecutable
sub skip_all_unless_exists {
    my ($type, $to_check) = @_;
    
    my $found;
    if ($type eq 'exe') {
        $found = can_execute($to_check);
    }
    elsif ($type eq 'env') {
        $found = $ENV{$to_check};
    }
    else {
        die "invalid type '$type'\n";
    }
    
    unless ($found) {
        my $skip_all = sub {
            my $builder = __PACKAGE__->builder;
            
            if ( not defined $builder->has_plan ) {
                $builder->skip_all(@_);
            }
            elsif ( $builder->has_plan eq 'no_plan' ) {
                $builder->skip(@_);
                if ( $builder->can('parent') && $builder->parent ) {
                    die bless {} => 'Test::Builder::Exception';
                }
                exit 0;
            }
            else {
                for ( 1 .. $builder->has_plan ) {
                    $builder->skip(@_);
                }
                if ( $builder->can('parent') && $builder->parent ) {
                    die bless {} => 'Test::Builder::Exception';
                }
                exit 0;
            }
        };
        
        if ($type eq 'exe') {
            $skip_all->("The test requires '$to_check' in PATH");
        }
        elsif ($type eq 'env') {
            $skip_all->("The test requires the '$to_check' environment variable");
        }
    }
}

sub can_execute {
    my $path = shift;
    
    if (is_file_path($path)) {
        return can_execute_path($path);
    }
    else {
        return which($path);
    }
}

sub can_execute_path {
    my $path = shift;
    if (-x $path) {
        if ($path =~ /\s/ && $path !~ /^$quote/) {
            $path = "$quote$path$quote";
        }
        return $path;
    }
    else {
        return;
    }
}

sub is_file_path {
    my $path = shift;
    return 1 if -e $path;
}

1;