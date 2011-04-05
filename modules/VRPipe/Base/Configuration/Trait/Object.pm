use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Object {
    use Data::Dumper;
    use VRPipe::Base::Configuration::Env;
    
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
    
    has _raw_config => (
        is      => 'ro',
        isa     => 'HashRef',
        lazy    => 1,
        builder => '_build_raw_config'
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
    
    method _build_raw_config {
        my $config_module = $self->config_module;
        try { eval "require $config_module"; die $@ if $@; } # require $config_module; with no eval does not work
        catch { return {} }
        return $config_module->get_config() || {};
    }
    
    method _from_config_or_env (Str $name, Str $env_key) {
        my $hash = $self->_raw_config();
        my $value = $hash->{$name};
        return $value if (defined $value && (ref($value) || $value ne '')); # $value may be an Env object which stringifies to '', which is allowed
        
        my $env;
        if (defined $ENV{lc $env_key}) {
            $env = lc $env_key;
        }
        elsif (defined $ENV{uc $env_key}) {
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
        
        my $values = {};
        my @options = $self->get_options;
        foreach my $option (@options) {
            my $key = $option->{key};
            my $value = $self->$key();
            next unless defined $value;
            $values->{$key} = $value;
        }
        
        my $dd = Data::Dumper->new( [ $values ], [ 'site_config' ] );
        open( my $fh, '>', $config_module_path ) or $self->throw("Cannot write to '$config_module_path': $!");
        printf $fh <<'END_HERE', $self->config_module, $dd->Dump();
package %s;
use strict;
use warnings;

my $site_config;
%s;

sub get_config {
    return $site_config;
}

1;
END_HERE
        close($fh);
    }
    
    method get_options ($class:) {
        my $meta = $class->meta;
        my @options;
        
        my @ck_attrs = sort { $a->question_number <=> $b->question_number } grep { $_->does('VRPipe::Base::Configuration::Trait::Attribute::ConfigKey') } $meta->get_all_attributes;
        foreach my $attr (@ck_attrs) {
            my $key = $attr->name;
            push(@options, { key => $key,
                             question => $attr->has_question ? $attr->question : "$key?",
                             question_number => $attr->question_number,
                             $attr->has_valid ? (valid => $attr->valid) : () });
        }
        return @options;
    }
}

1;