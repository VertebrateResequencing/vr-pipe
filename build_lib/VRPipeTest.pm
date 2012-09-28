
=head1 NAME

VRPipeTest - for use by all test scripts

=head1 SYNOPSIS
    
    BEGIN {
        use Test::Most tests => 3;
        use VRPipeTest (required_env => [qw(RUN_AUTHOR_TESTS)],
                        required_exe => [qw(foo bar)],
                        required_module => [qw(Baz::Bar Far::Car)],
                        debug => 0,
                        max_retries => 3);
    }

Set debug to 1 and max_retries to 0 when you're trying to investigate really
obscure issues with a pipeline.

=head1 DESCRIPTION

*** more documentation to come

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

package VRPipeTest;
use strict;
use warnings;
use base 'Test::DBIx::Class';
use VRPipe::Persistent::SchemaBase;
use File::Spec;
use File::Which qw(which);
$SQL::Translator::Schema::DEBUG = 0; # suppress stupid warning in test harness
$ENV{PERL_STRICTURES_EXTRA} = 0;

BEGIN { unshift(@INC, './modules', './t') }
VRPipe::Persistent::SchemaBase->database_deployment('testing'); # do not change this under any circumstance!

use constant WIN32 => $^O eq 'MSWin32';
my $quote = WIN32 ? q/"/ : q/'/;
our $max_retries = 3;
our $debug       = 0;

sub import {
    my $class  = shift;
    my $caller = caller(0);
    my %args   = @_;
    
    # skip if required_env or required_exe not satisfied
    while (my ($key, $val) = each %args) {
        if ($key =~ /^required_(env|exe|module)/) {
            my $type = $1;
            my @to_check = ref($val) ? @{$val} : ($val);
            foreach (@to_check) {
                skip_all_unless_exists($type => $_);
            }
        }
        elsif ($key eq 'max_retries') {
            $max_retries = $val;
        }
        elsif ($key eq 'debug') {
            $debug = $val;
        }
    }
    
    # copy/pasted  Test::DBIx::Class to get the exports
    my ($schema_manager, $merged_config, @exports) = $class->_initialize();
    my $exporter = Sub::Exporter::build_exporter({
            exports => [
                Schema => sub {
                    return sub {
                        return $schema_manager->schema;
                    };
                },
                get_elements => sub {
                    return \&get_elements;
                },
                result_with_inflated_paths => sub {
                    return \&result_with_inflated_paths;
                  }
            ],
            into_level => 1,
        }
    );
    
    $class->$exporter(qw/Schema get_elements result_with_inflated_paths/, @exports);
}

sub _initialize_schema {
    my $class = shift;
    
    my $config = {
        schema_class     => 'VRPipe::Persistent::Schema',
        force_drop_table => 1,
        keep_db          => 1,
        connect_info     => [VRPipe::Persistent::SchemaBase->get_dsn, VRPipe::Persistent::SchemaBase->get_user, VRPipe::Persistent::SchemaBase->get_password]
    };
    
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
    elsif ($type eq 'module') {
        eval "use $to_check;";
        $found = !$@;
    }
    else {
        die "invalid type '$type'\n";
    }
    
    unless ($found) {
        my $skip_all = sub {
            my $builder = __PACKAGE__->builder;
            
            if (not defined $builder->has_plan) {
                $builder->skip_all(@_);
            }
            elsif ($builder->has_plan eq 'no_plan') {
                $builder->skip(@_);
                if ($builder->can('parent') && $builder->parent) {
                    die bless {} => 'Test::Builder::Exception';
                }
                exit 0;
            }
            else {
                for (1 .. $builder->has_plan) {
                    $builder->skip(@_);
                }
                if ($builder->can('parent') && $builder->parent) {
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
        elsif ($type eq 'module') {
            $skip_all->("The test requires the '$to_check' Perl module");
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

sub debug {
    return $debug;
}

sub max_retries {
    return $max_retries;
}

sub get_elements {
    my $ds       = shift;
    my $pager    = $ds->elements;
    my @elements = ();
    while (my $elements = $pager->next) {
        push(@elements, @$elements);
    }
    return \@elements;
}

sub result_with_inflated_paths {
    my $de     = shift;
    my $result = $de->result;
    if (defined $result->{paths}) {
        $result->{paths} = [$de->paths] unless ref($result->{paths});
    }
    return $result;
}

1;
