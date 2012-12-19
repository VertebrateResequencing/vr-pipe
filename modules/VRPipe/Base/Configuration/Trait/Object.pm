
=head1 NAME

VRPipe::Base::Configuration::Trait::Object - config object methods

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

role VRPipe::Base::Configuration::Trait::Object {
    use Data::Dumper;
    use VRPipe::Base::Configuration::Env;
    use VRPipe::Base::Configuration::Option;
    use Crypt::CBC;
    use Digest::MD5 qw(md5_hex);
    
    has config_module => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_config_module'
    );
    
    has config_module_path => (
        is      => 'ro',
        isa     => File,
        coerce  => 1,
        lazy    => 1,
        builder => '_build_config_path'
    );
    
    has _key_file => (
        is      => 'rw',
        isa     => MaybeFile,
        coerce  => 1,
        lazy    => 1,
        builder => '_build_key_file'
    );
    
    has _raw_config => (
        is      => 'ro',
        isa     => 'HashRef',
        lazy    => 1,
        builder => '_build_raw_config'
    );
    
    has _options_list => (
        is      => 'ro',
        isa     => 'ArrayRef',
        builder => '_build_options_list'
    );
    
    has _next_option_index => (
        is      => 'rw',
        isa     => 'Int',
        default => 0
    );
    
    method _build_config_module {
        return 'VRPipe::SiteConfig';
    }
    
    method _build_config_path {
        my $module = $self->config_module;
        
        my $path = file(split('::', $module));
        $path .= '.pm';
        
        eval "require $module";
        if (exists $INC{$path}) {
            $path = $INC{$path};
        }
        else {
            # stick it in the first writable path in INC
            my $writable_inc;
            foreach my $inc (@INC) {
                if (-w $inc) { # only a permission test; doesn't really guarantee we can really write there
                    $writable_inc = $inc;
                    last;
                }
            }
            $writable_inc || $self->throw("There don't seem to be any writable paths in your INC");
            
            $path = file($writable_inc, $path);
        }
        
        return $path;
    }
    
    method _build_key_file {
        if ($self->can('encryption_key_file')) {
            return $self->encryption_key_file;
        }
        return;
    }
    
    method _build_raw_config {
        my $config_module = $self->config_module;
        try { eval "require $config_module"; die $@ if $@; } # require $config_module; with no eval does not work
        catch { return {} } my $config = $config_module->get_config() || {};
        
        # decrypt any encrypted values
        my $key_file;
        if (exists $config->{encryption_key_file}) {
            $key_file = $self->_key_file($config->{encryption_key_file});
        }
        while (my ($key, $value) = each %{$config}) {
            if ($value =~ /^Salted__/) {                     #*** How stable is this 'Salted__' string? Should probably check to see if the corresponding attribute is secure instead
                if ($key_file && -s $key_file) {
                    my $fh             = $key_file->openr();
                    my $decryption_key = <$fh>;
                    $fh->close;
                    chomp($decryption_key);
                    
                    my $cipher = Crypt::CBC->new(
                        -key    => $decryption_key,
                        -cipher => 'Blowfish',
                        -header => 'salt'
                    );
                    
                    $config->{$key} = $cipher->decrypt($value);
                }
                else {
                    $self->throw("Cannot read an encrypted config option ($key) when the encryption_key_file option is not set to a non-empty file");
                }
            }
        }
        
        return $config;
    }
    
    method _from_config_or_env (Str $name, Str $env_key) {
        my $hash  = $self->_raw_config();
        my $value = $hash->{$name};
        return $value if (defined $value && (ref($value) || $value ne '')); # $value may be an Env object which stringifies to '', which is allowed
        
        my $env;
        if (defined $ENV{ lc $env_key }) {
            $env = lc $env_key;
        }
        elsif (defined $ENV{ uc $env_key }) {
            $env = uc $env_key;
        }
        if ($env) {
            # return an env object that stringifies to $ENV{$env}, but lets us
            # store the fact that we're dealing with an environment variable
            # in the config file
            return VRPipe::Base::Configuration::Env->new(variable => $env_key);
        }
        
        return;
    }
    
    method write_config_module {
        my $config_module_path = $self->config_module_path();
        
        my $values            = {};
        my $next_option_index = $self->_next_option_index;
        $self->_next_option_index(0);
        while (my $option = $self->next_option) {
            my $value = $option->value();
            next unless defined $value;
            
            if ($option->secure && !ref($value)) {
                my $key_file = $self->_key_file;
                unless ($key_file) {
                    $self->throw("Secure config options cannot be set without the encryption_key_file option being set");
                }
                
                my $key;
                if (-s $key_file) {
                    my $fh = $key_file->openr();
                    $key = <$fh>;
                    $fh->close;
                    unless (defined $key) {
                        $self->throw("The encryption key file you specified ($key_file) has something in it, but isn't readable; please make it readable by yourself and the user that will run vrpipe-server.");
                    }
                    unless (length($key) == 33) {
                        $self->throw("The encryption key file you specified ($key_file) has something in it, but the first line doesn't seem to specify a valid key - is it really a key file?");
                    }
                    chomp($key);
                }
                else {
                    $key = md5_hex(rand);
                    my $perms = sprintf "%o", "600";
                    my $fh = $key_file->open('w', $perms);
                    unless ($fh) {
                        $self->throw("Could not open '$key_file': $!; $@");
                    }
                    print $fh $key, "\n";
                    $fh->close;
                }
                
                my $cipher = Crypt::CBC->new(
                    -key    => $key,
                    -cipher => 'Blowfish',
                    -header => 'salt'
                );
                
                $value = $cipher->encrypt($value);
            }
            
            $values->{ $option->key } = $value;
        }
        $self->_next_option_index($next_option_index);
        
        my $dd = Data::Dumper->new([$values], ['config']);
        open(my $fh, '>', $config_module_path) or $self->throw("Cannot write to '$config_module_path': $!");
        printf $fh <<'END_HERE', $self->config_module, $dd->Dump();
package %s;
use strict;
use warnings;

my $config;
%s

sub get_config {
    return $config;
}

1;
END_HERE
        close($fh);
    }
    
    method _build_options_list ($class:) {
        my $meta = $class->meta;
        my @ck_attrs = sort { $a->question_number <=> $b->question_number } grep { $_->does('VRPipe::Base::Configuration::Trait::Attribute::ConfigKey') } $meta->get_all_attributes;
        
        my $do_add = 0;
        foreach my $attr (@ck_attrs) {
            my $name = $attr->name;
            if ($name eq 'encryption_key_file') {
                $do_add = 0;
                last;
            }
            if ($attr->has_secure) {
                $do_add = 1;
            }
        }
        
        if ($do_add) {
            # add a question about where the encryption key should be stored
            $meta->make_mutable;
            my $new_attr = $meta->add_attribute(
                'encryption_key_file',
                accessor        => 'encryption_key_file',
                is              => 'rw',
                trigger         => sub { shift->_key_file(shift); },
                question        => 'Passwords you enter will be encrypted; where should your encryption key be stored? (it is up to you to properly secure this file)',
                question_number => 0
            );
            $meta->make_immutable;
            unshift(@ck_attrs, $new_attr);
        }
        
        return \@ck_attrs;
    }
    
    method next_option {
        my $next_option_index = $self->_next_option_index;
        my $attrs             = $self->_options_list;
        my $max               = $#$attrs;
        if ($next_option_index > $max) {
            return;
        }
        
        my $attr;
        while (!defined $attr) {
            $attr = $attrs->[$next_option_index++];
            if ($attr->has_skip) {
                my $skip_method = $attr->skip;
                if ($self->$skip_method) {
                    undef $attr;
                }
            }
            last if $next_option_index > $max;
        }
        
        $self->_next_option_index($next_option_index);
        unless ($attr) {
            return;
        }
        return $self->_new_option($attr);
    }
    
    method option (Int :$number?, Str :$key?) {
        my $attrs = $self->_options_list;
        
        my $attr;
        if ($number) {
            foreach my $this_attr (@$attrs) {
                if ($this_attr->question_number == $number) {
                    $attr = $this_attr;
                    last;
                }
            }
        }
        elsif ($key) {
            foreach my $this_attr (@$attrs) {
                if ($this_attr->name eq $key) {
                    $attr = $this_attr;
                    last;
                }
            }
        }
        
        unless (defined $attr) {
            $self->throw("You requested an option that doesn't exist");
        }
        
        return $self->_new_option($attr);
    }
    
    method _new_option ($attr) {
        return VRPipe::Base::Configuration::Option->new(_attr => $attr, _obj => $self);
    }
}

1;
